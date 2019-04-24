/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <cmath>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1.
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method
   *   (and others in this file).
   */
  num_particles = 1000;  // TODO: Set the number of particles

  is_initialized = true;

  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(fmod(theta, (2*M_PI)), std[2]);

  for (int i=0; i<num_particles; ++i)
  {
    weights.push_back(1.0);

    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    particles.push_back(p);
  }

}

void ParticleFilter::prediction(double delta_t, double std_pos[],
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

   std::default_random_engine gen;
   // normal_distribution<double> dist_x(0.0, std[0]);
   // normal_distribution<double> dist_y(0.0, std[1]);
   // normal_distribution<double> dist_theta(0.0, std[2]);

   for (int i=0; i<num_particles; ++i)
   {
     double x_0 = particles[i].x;
     double y_0 = particles[i].y;
     double theta_0 = fmod(particles[i].theta, (2*M_PI));

     double theta_f = fmod((theta_0 + yaw_rate*delta_t), (2*M_PI));

     // double new_theta = theta_f + dist_theta(gen);
     normal_distribution<double> dist_theta(theta_f, std_pos[2]);
     double new_theta = dist_theta(gen);
     particles[i].theta = new_theta;

     double dx = (velocity/yaw_rate)*(sin(theta_f) - sin(theta_0));
     double x_f = x_0 + dx;
     double dy = (velocity/yaw_rate)*(cos(theta_0) - cos(theta_f));
     double y_f = y_0 + dy;

     // double new_x = x_f + dist_x(gen);
     normal_distribution<double> dist_x(x_f, std_pos[0]);
     double new_x = dist_x(gen);
     particles[i].x = new_x;

     // double new_y = y_f + dist_y(gen);
     normal_distribution<double> dist_y(y_f, std_pos[1]);
     double new_y = dist_y(gen);
     particles[i].y = new_y;

   }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted,
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each
   *   observed measurement and assign the observed measurement to this
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will
   *   probably find it useful to implement this method and use it as a helper
   *   during the updateWeights phase.
   */

   for (int i=0; i<observations.size(); ++i)
   {
     if (!predicted.empty())
     {
       double min_d = dist(predicted[0].x, predicted[0].y,
                            observations[i].x, observations[i].y);
       observations[i].id = predicted[0].id;
       int min_i = 0;

       for (int j=1; j<predicted.size(); ++j)
       {
         double d = dist(predicted[j].x, predicted[j].y,
                          observations[i].x, observations[i].y);

         if (d < min_d)
         {
           observations[i].id = predicted[j].id;
           min_d = d;
           min_i = j;
         }
       }

       predicted.erase(predicted.begin()+min_i);

     }

   }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const vector<LandmarkObs> &observations,
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian
   *   distribution. You can read more about this distribution here:
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system.
   *   Your particles are located according to the MAP'S coordinate system.
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

   for (int i=0; i<num_particles; ++i)
   {

     double x_part = particles[i].x;
     double y_part = particles[i].y;
     double theta = fmod(particles[i].theta, (2*M_PI));

     vector<LandmarkObs> observations_trans;
     for (int j=0; j<observations.size(); ++j)
     {
       double x_obs = observations[j].x;
       double y_obs = observations[j].y;

       LandmarkObs observation;
       observation.x = x_part + (cos(theta)*x_obs) - (sin(theta)*y_obs);
       observation.y = y_part + (sin(theta)*x_obs) + (cos(theta)*y_obs);

       observations_trans.push_back(observation);
     }

     vector<LandmarkObs> predicted;
     for (int j=0; j<map_landmarks.landmark_list.size(); ++j)
     {
       double x_landmark = map_landmarks.landmark_list[j].x_f;
       double y_landmark = map_landmarks.landmark_list[j].y_f;
       double d = dist(x_landmark, y_landmark, x_part, y_part);
       if (d<sensor_range)
       {
         LandmarkObs landmark;
         landmark.id = map_landmarks.landmark_list[j].id_i;
         landmark.x = x_landmark;
         landmark.y = y_landmark;

         predicted.push_back(landmark);
       }
     }

     dataAssociation(predicted, observations_trans);

     vector<int> associations;
     vector<double> sense_x;
     vector<double> sense_y;

     double tot_prob = 1.0;

     double sigma_x = std_landmark[0];
     double sigma_y = std_landmark[1];
     double gauss_norm = 1/(2*M_PI*sigma_x*sigma_y);
     for (int j=0; j<observations_trans.size(); ++j)
     {

       double x_obs = observations_trans[j].x;
       double y_obs = observations_trans[j].y;
       int observe_id = observations_trans[j].id;

       for (int k=0; k<predicted.size(); ++k)
       {
         int predicted_id = predicted[k].id;
         if (predicted_id==observe_id)
         {

           double mu_x = predicted[k].x;
           double mu_y = predicted[k].y;

           associations.push_back(predicted_id);
           sense_x.push_back(mu_x);
           sense_y.push_back(mu_y);

           double exponent = (pow(x_obs - mu_x, 2)/(2*pow(sigma_x, 2)))
                              + (pow(y_obs-mu_y, 2)/(2*pow(sigma_y, 2)));

           double gauss_prob = gauss_norm*exp(-exponent);

           tot_prob *= gauss_prob;

         }
       }
     }

     particles[i].weight = tot_prob;
     weights[i] = tot_prob;

     SetAssociations(particles[i], associations, sense_x, sense_y);

     observations_trans.clear();
     predicted.clear();
   }

   double weight_sum = 0;
   for (int i=0; i<num_particles; ++i)
   {
     weight_sum += weights[i];
   }

   for (int i=0; i<num_particles; ++i)
   {
     particles[i].weight /= weight_sum;
     weights[i] /= weight_sum;
   }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional
   *   to their weight.
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

   std::default_random_engine gen;
   std::discrete_distribution<> resample_dist(weights.begin(), weights.end());

   vector<Particle> particles_resample;
   for (int i=0; i<num_particles; ++i)
   {
     int rand_i = resample_dist(gen);
     particles_resample.push_back(particles[rand_i]);
   }
   particles = particles_resample;

}

void ParticleFilter::SetAssociations(Particle& particle,
                                     const vector<int>& associations,
                                     const vector<double>& sense_x,
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association,
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
