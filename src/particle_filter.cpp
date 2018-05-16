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
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	// Number of particles to draw:
	num_particles = 100;

	// Random engine init:
	default_random_engine gen;

	// Setting std to x, y, theta:
	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];

	// Creating normal distributions:
	std::normal_distribution<double> dist_x(x, std_x);
	std::normal_distribution<double> dist_y(y, std_y);
	std::normal_distribution<double> dist_theta(theta, std_theta);

	// Looping through the number of particles and generating randomnes:
	for (int i=0; i<num_particles; i++) {
		Particle P;
		P.id = i;
		P.x = dist_x(gen);
		P.y = dist_y(gen);
		P.theta = dist_theta(gen);
		P.weight = 1;

		particles.push_back(P);
		weights.push_back(P.weight);
	}
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	// Random engine init:
	default_random_engine gen;
	// Initializing new variables:
	double new_x;
	double new_y;
	double new_theta;

	// Looping through the particles:
	for (int i=0; i<num_particles; i++) {
		// Checking whether or not yaw_rate is 0:
		if (fabs(yaw_rate) < 0.0001) {
			// Updating variables if yaw_rate is 0:
			new_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			new_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
			new_theta = particles[i].theta;
		} else {
			// Updating variables if yaw_rate is NOT 0:
			new_x = particles[i].x + (velocity / yaw_rate) *
					(sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
			new_y = particles[i].y + (velocity / yaw_rate) *
					(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			new_theta = particles[i].theta + yaw_rate * delta_t;
		}

		// Creating normal distributions:
		normal_distribution<double> dist_x(new_x, std_pos[0]);
		normal_distribution<double> dist_y(new_y, std_pos[1]);
		normal_distribution<double> dist_theta(new_theta, std_pos[2]);

		// Particles and generating randomnes
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.
	// Looping through the observations:
	for (unsigned int i=0; i<observations.size(); i++) {
		// Variable that will keep track of the closest landmarks:
		int closest_landmark = -1;
		// Variable that keeps track of the minimun distance:
		double min_distance = 100000;

		// Looping through the predictions:
		for (unsigned int j=0; j<predicted.size(); j++) {
			// Variable that collects the observations and predictions of x & y:
			double distance = dist(observations[i].x,
									observations[i].y,
									predicted[j].x,
									predicted[j].y);
			// Checking whether or not the distance is lower than the min_distance previously set:
			if (distance < min_distance) {
				// Updating the min_distance if the distance is smaller;
				min_distance = distance;
				// Saving the object to a variable:
				closest_landmark = predicted[j].id;
			}
		}
		// Updating the observations whether or not the distances change:
		observations[i].id = closest_landmark;
	}
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
	//   http://planning.cs.uiuc.edu/node99.html
	double weight_normalizer = 0.0;

	for (int i=0; i<num_particles; i++) {
	    // Creating a vector to hold the map landmark locations predicted to be within sensor range:
	    vector<LandmarkObs> filtered_map;
		double distance;

	    // Looping through the landmarks:
	    for (unsigned int j=0; j<map_landmarks.landmark_list.size(); j++) {
			// Saving id and x, y coordinates:
			distance = dist(particles[i].x,
							particles[i].y,
							map_landmarks.landmark_list[j].x_f,
							map_landmarks.landmark_list[j].y_f);

			// only consider landmarks within sensor range of the particle (rather than using the "dist" method considering a circular
			if (distance < sensor_range) {
				// Adding prediction to vector:
				LandmarkObs obs;
				obs.id = map_landmarks.landmark_list[j].id_i;
				obs.x = map_landmarks.landmark_list[j].x_f;
				obs.y = map_landmarks.landmark_list[j].y_f;
				filtered_map.push_back(obs);
			}
	    }

	    // Creating and populating a copy of the list of observations transformed from vehicle coordinates:
	    vector<LandmarkObs> transformed_obs;
		// Looping through observations:
	    for (unsigned int j=0; j<observations.size(); j++) {
			LandmarkObs obs;
		    obs.id = j;
		    obs.x = particles[i].x + (cos(particles[i].theta) * observations[j].x) -
									(sin(particles[i].theta) * observations[j].y);
		    obs.y = particles[i].y + (sin(particles[i].theta) * observations[j].x) +
			 						(cos(particles[i].theta) * observations[j].y);
		    transformed_obs.push_back(obs);
	    }

	    // Calling dataAssociation function on filtered_map and transformed_os
	    dataAssociation(filtered_map, transformed_obs);

	    // Re-initializing weight:
	    particles[i].weight = 1.0;
		double sigma_x = std_landmark[0];
		double sigma_y = std_landmark[1];
		double sigma_x_2 = pow(sigma_x, 2);
		double sigma_y_2 = pow(sigma_y, 2);
		double normalizer = (1.0/(2.0 * M_PI * sigma_x * sigma_y));
		int k, l;

		// Looping through the transformed_os:
	    for (k=0; k<transformed_obs.size(); k++) {
			// placeholders for observation and associated prediction coordinates
			double trans_obs_x = transformed_obs[k].x;
			double trans_obs_y = transformed_obs[k].y;
			double trans_obs_id = transformed_obs[k].id;
			double multi_prob = 1.0;

			// get the x,y coordinates of the prediction associated with the current observation
			for (l=0; l<filtered_map.size(); l++) {
				double pred_landmark_x = filtered_map[l].x;
				double pred_landmark_y = filtered_map[l].y;
				double pred_landmark_id = filtered_map[l].id;

				if (trans_obs_id == pred_landmark_id) {
					// Calculating weights for observations with multivariate Gaussian:
					multi_prob = normalizer * exp(-1.0 * (
						(pow((trans_obs_x - pred_landmark_x), 2)/(2.0 * sigma_x_2)) + (pow((trans_obs_y - pred_landmark_y), 2)/(2.0 * sigma_y_2))));
					particles[i].weight *= multi_prob;
				}
			}
	    }
		weight_normalizer += particles[i].weight;
	}
	// Looping through the updated weights:
	for (int i = 0; i < particles.size(); i++) {
		particles[i].weight /= weight_normalizer;
		weights[i] = particles[i].weight;
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// Creating new vector for resampled particles:
	vector<Particle> resampled_particles;
	// Random engine init:
	default_random_engine gen;
	// Random particle index:
	uniform_int_distribution<int> particle_index(0, num_particles - 1);
	// Generating random particle index:
	int current_index = particle_index(gen);
	// Initializing a variable at 0 to keep track of weights:
	double beta = 0.0;
	// Setting the maximun weights of the particles:
	double max_weight = 2.0 * *max_element(weights.begin(), weights.end());

	// Looping through the particles:
	for (int i=0; i<particles.size(); i++) {
		// Creating variable and setting parameters:
		uniform_real_distribution<double> random_weight(0.0, max_weight);
		// Giving all the particles randome weights:
		beta += random_weight(gen);

		// Looping until the beta is smaller than the weights:
		while (beta > weights[current_index]) {
			// Updating the beta result:
			beta -= weights[current_index];
			// Updating the index of the weights:
			current_index = (current_index + 1) % num_particles;
		}
	// Resampling the data:
	resampled_particles.push_back(particles[current_index]);
	}
	// Updating the particles:
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
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
