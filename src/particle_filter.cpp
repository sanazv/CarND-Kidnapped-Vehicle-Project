/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    num_particles_ = 10;
    particles_.reserve(num_particles_);
    
    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);
    
    for (int i=0; i<num_particles_; ++i)
    {
        Particle p(i+1, dist_x(gen), dist_y(gen), dist_theta(gen));
        particles_.push_back(p);
    }
    
    is_initialized_ = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    
    default_random_engine gen;
    normal_distribution<double> dist_x(0.0, sqrt(delta_t)*std_pos[0]);
    normal_distribution<double> dist_y(0.0, sqrt(delta_t)*std_pos[1]);
    normal_distribution<double> dist_theta(0.0, sqrt(delta_t)*std_pos[2]);
    
    for (int i=0; i<num_particles_; ++i)
    {
        if (abs(yaw_rate)>0.00001)
        {
            particles_[i].x += (velocity/yaw_rate)*(sin(particles_[i].theta + yaw_rate*delta_t) - sin(particles_[i].theta)) + dist_x(gen);
            particles_[i].y += (velocity/yaw_rate)*(cos(particles_[i].theta) - cos(particles_[i].theta + yaw_rate*delta_t) ) + dist_y(gen);
            particles_[i].theta += yaw_rate*delta_t + dist_theta(gen);
        } else {
            particles_[i].x +=  velocity*cos(particles_[i].theta)*delta_t + dist_x(gen);
            particles_[i].y +=  velocity*sin(particles_[i].theta)*delta_t + dist_y(gen);
            particles_[i].theta += dist_theta(gen);
        }

    }
    

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    double weight_sum = 0.0;
    
    // for each particle, let's go through the observations and map them to the landmarks
    for (int i=0; i<particles_.size(); ++i)
    {
        Particle & p = particles_[i];
        p.associations.clear();
        p.sense_x.clear();
        p.sense_y.clear();
        
        double prob = 1.0;
        
        // Given the noisy observations, we need to translate them to the map coords
        for (int obs = 0; obs<observations.size(); ++obs)
        {
            LandmarkObs & lo = observations[obs];
            
            // TODO convert the observation to the map coordinates
            double x_map = p.x + lo.x*cos(p.theta) - lo.y*sin(p.theta);
            double y_map = p.y + lo.x*sin(p.theta) + lo.y*cos(p.theta);
            
            // calculate the distance to each of the known landmarks, also find the nearest
            double minDistancesToLandmarks = 1.0e10;
            int minDistancesToLandmarksID = -1;
            for (int l = 0; l<map_landmarks.landmark_list.size(); ++l)
            {
                // note that dist is an expensive function, we can use the sensor_range to reduce the number of calls
                double d = dist(x_map, y_map, map_landmarks.landmark_list[l].x_f, map_landmarks.landmark_list[l].y_f);
                if (d<minDistancesToLandmarks)
                {
                    minDistancesToLandmarks = d;
                    minDistancesToLandmarksID = l;
                }
            }
            
            // compute the weight as the liklihood of the current particle position given the measurements.
            double denom = 1.0/(2.0*M_PI*std_landmark[0]*std_landmark[1]);
            double t1 = pow(x_map - map_landmarks.landmark_list[minDistancesToLandmarksID].x_f,2) / (std_landmark[0]*std_landmark[0]);
            double t2 = pow(y_map - map_landmarks.landmark_list[minDistancesToLandmarksID].y_f,2) / (std_landmark[1]*std_landmark[1]);
            prob *= denom * exp(-0.5*( t1 + t2 ));

        }
        
        p.weight = prob;
        weight_sum += prob;
    }
    
    // print the weights
    //for (int i=0; i<particles_.size(); ++i)
    //    cout << i << "," << particles_[i].weight << endl;
    
    // normalize the weights vector
    double sum = 0.0;
    for (int i=0; i<particles_.size(); ++i)
        sum += particles_[i].weight;
    
    cout << "sum: " << sum << endl;
    
    for (int i=0; i<particles_.size(); ++i)
        particles_[i].weight /= sum;
    
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    
    int N = particles_.size();
    vector<double> weights;
    for (int i=0; i<particles_.size(); ++i)
        weights.push_back(particles_[i].weight);
    
    double w_max = -1.0;
    for (int i = 0; i < N; ++i) {
        if (weights[i] > w_max) {
            w_max = weights[i];
        }
    }
    
    
    double beta = 0.0;
    std::vector<Particle> temp_particles;
    
    int index = rand()%N;
    
    std::default_random_engine generator;
    std::uniform_real_distribution<double> u_distrib(0.0,1.0);
    
    for (int i = 0; i < N; ++i) {
        beta = beta + 2 * w_max * u_distrib(generator);
        
        while (weights[index] < beta){
            beta = beta - weights[index];
            index = index + 1;
            index = index%N;
        }
        temp_particles.push_back(particles_[index]);
    }
    
    particles_ = temp_particles;
    
    
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
