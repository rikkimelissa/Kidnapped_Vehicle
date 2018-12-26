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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult partnicle_filter.h for more information about this method (and others in this file).
	num_particles = 5;  // TODO: Set the number of particles

    default_random_engine gen;
    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    for (int i=0; i<num_particles; i++){
        Particle newparticle = Particle();
        newparticle.id = i;
        newparticle.x = dist_x(gen);
        newparticle.y = dist_y(gen);
        newparticle.theta = dist_theta(gen);
        newparticle.weight = 1;
        particles.push_back(newparticle);
    }

    is_initialized = 1;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;

    // std::cout << << endl;
    // std::cout << velocity << endl;
    // std::cout << yaw_rate << endl;


    for (int i=0; i<num_particles; i++){
        normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
        normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
        normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

        particles[i].x = dist_x(gen) + velocity/yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
        particles[i].y = velocity/yaw_rate*(-cos(particles[i].theta + yaw_rate*delta_t) + cos(particles[i].theta)) + dist_y(gen);
        particles[i].theta = yaw_rate*delta_t + dist_theta(gen);
    }

    // std::cout << "finished predicting" << endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    double distMin, dist;

    for (int o=0; o<observations.size(); o++){
        distMin = 1000000;
        for (int p=0; p<predicted.size(); p++){
            dist = sqrt((predicted[p].x-observations[o].x)*(predicted[p].x-observations[o].x)+(predicted[p].y-observations[o].y)*(predicted[p].y-observations[o].y));
            if (dist < distMin){
                distMin = dist;
                observations[o].id = predicted[p].id;
            }
        }    
    }

    // std::cout << "finished data association" << endl;

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.htmlupdateWeigh

    double ux, uy;
    double x, y, theta;
    std::vector<double> weights;
    LandmarkObs newLandmark;

    for (int p=0; p<num_particles; p++){ 

        x = particles[p].x;
        y = particles[p].y;
        theta = particles[p].theta;
        std::vector<LandmarkObs> vis_landmarks;
        std::vector<LandmarkObs> obs_landmarks;

        std::cout << p << std::endl;
        std::cout << x << std::endl;
        std::cout << y << std::endl;

        // Find landmarks in visible range of particle
        for (int l=0; l<map_landmarks.landmark_list.size(); l++){
            if (sqrt((map_landmarks.landmark_list[l].x_f - x)*(map_landmarks.landmark_list[l].x_f - x) + (map_landmarks.landmark_list[l].y_f - y)*(map_landmarks.landmark_list[l].y_f - y)) < sensor_range){
                newLandmark.x = map_landmarks.landmark_list[l].x_f;
                newLandmark.y = map_landmarks.landmark_list[l].y_f;
                newLandmark.id = map_landmarks.landmark_list[l].id_i;
                vis_landmarks.push_back(newLandmark);
            }
        }

        // std::cout << "finished finding landmarks" << endl;


        // Transform observations to map coordinates
        for (int i=0; i<observations.size(); i++){
            newLandmark.x = x + (cos(theta)*observations[i].x - sin(theta)*observations[i].y);
            newLandmark.y = y + (sin(theta)*observations[i].x + cos(theta)*observations[i].y);
            obs_landmarks.push_back(newLandmark);
        }


        // Associate observed landmarks with map landmark
        dataAssociation(vis_landmarks, obs_landmarks);

        // Update particle weights using multi-variate Gaussian probability density function
        double prob = 1.0;
        double gauss_norm, exponent, weight;
        for (int l=0; l<obs_landmarks.size(); l++){
            ux = map_landmarks.landmark_list[obs_landmarks[l].id-1].x_f;
            uy = map_landmarks.landmark_list[obs_landmarks[l].id-1].y_f;
            gauss_norm = 1.0/(2*acos(-1.0)*std_landmark[0]*std_landmark[1]);
            exponent = ((obs_landmarks[l].x-ux)*(obs_landmarks[l].x-ux)/(2*std_landmark[0]*std_landmark[0])+(obs_landmarks[l].y-uy)*(obs_landmarks[l].y-uy)/(2*std_landmark[1]*std_landmark[1]));
            weight = gauss_norm*exp(-exponent);
            prob *= weight;
        }

        particles[p].weight = prob;
        weights.push_back(prob);
    }

    std::cout << "normalizing" << endl;
    double maxW = *max_element(weights.begin(), weights.end());
    std::cout << maxW << endl;

    // Normalize weights from 0 to 1
    for (int p=0; p<num_particles; p++){
        // std::cout << particles[p].x << std::endl;
        // std::cout << particles[p].y << std::endl;
        std::cout << particles[p].weight << std::endl; 
        particles[p].weight /= maxW;
    }

    // std::cout << "finished normalizing" << endl;

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    std::vector<Particle> updated_particles;
    int index = 0;
    float beta = 0;

    float randNum;

    // std::cout << "resampled" << std::endl;

    for (int i=0; i<num_particles; i++){
        randNum = rand()/double(RAND_MAX);
        beta = beta + randNum;

        while (particles[index].weight < beta){
            beta = beta - particles[index].weight;
            index += 1;
            if (index > (num_particles-1)){
                index = 0;
            }
        }
        updated_particles.push_back(particles[index]);
        std::cout << index << endl;
    }

    // std::cout << "finished resampling" << endl;
    particles = updated_particles;

}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
