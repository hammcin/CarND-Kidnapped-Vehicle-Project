# Kidnapped Vehicle Project

### Goals

In this project you will implement a 2 dimensional particle filter in C++. Your particle filter will be given a map and some initial localization information (analogous to what a GPS would provide). At each time step your filter will also get observation and control data.

### Accuracy and Performance

#### 1. Does your particle filter localize the vehicle to within the desired accuracy and run within the specified time of 100 seconds?

Here's a [link to my video result](./kidnapped_vehicle_video.mp4).

### General

#### 1. Does your code use a particle filter to localize the robot?

A particle filter is used to localize the kidnapped vehicle (particle_filter.h and particle_filter.cpp).  The coordinates of the particles in the filter are initialized using a GPS estimate with Gaussian noise (particle_filter.cpp, ParticleFilter::init, lines 26-57).

On each iteration of the filter, the next state of the vehicle is predicted using the velocity and yaw rate from the previous iteration (particle_filter.cpp, ParticleFilter::prediction, lines 59-114).  Prediction is performed using the bicycle model for the case of zero yaw rate (particle_filter.cpp, ParticleFilter::prediction, lines 81-88) and non-zero yaw rate (particle_filter.cpp, ParticleFilter::prediction, lines 89-98).

The weights of each particle are then updated (particle_filter.cpp, ParticleFilter::updateWeights, lines 157-274).  The coordinates of the landmark observations must first be transformed from the vehicle's coordinate system to the map coordinate system (particle_filter.cpp, ParticleFilter::updateWeights, lines 177-192).  Associations between landmark measurements and landmark identification numbers are obtained using nearest neighbors (particle_filter.cpp, ParticleFilter::dataAssociation, lines 116-155).  The weights are then updated based on the likelihood of the landmark measurements for each particle using a multi-variate Gaussian distribution (particle_filter.cpp, ParticleFilter::updateWeights, lines 217-251).  Finally, the particle weights are normalized (particle_filter.cpp, ParticleFilter::updateWeights, lines 262-272).

The weight update step is followed by resampling from the set of particles (particle_filter.cpp, ParticleFilter::resample, lines 276-295).  Particles are resampled in such a way that the probability of sampling a particle is proportional to the particle's weight as computed in the weight update step.  This is accomplished by using the C++ library's `std::discrete_distribution`.
