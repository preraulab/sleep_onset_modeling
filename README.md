# Tracking The Sleep Onset Process (SOP)
#### In this repository, we develop new way of tracking the continuum of sleep onset, which incorporates EEG, physiological, and behavioral information to estimate the instantaneous probability that a subject is awake at each point in time.
## Citation
> Prerau, M. J., Hartnack, K. E., Obregon-Henao, G., Sampson, A., Merlino, M., Gannon, K., Bianchi, M., Ellenbogen, J., & Purdon, P. L. (2014). Tracking the sleep onset process: an empirical model of behavioral and physiological dynamics. PLoS computational biology, 10(10), e1003866. https://doi.org/10.1371/journal.pcbi.1003866
<br/><br/>
## Overview
<p align="center"> 
<img src="https://prerau.bwh.harvard.edu/wp-content/uploads/2022/10/SOP.png" alt="sleeop onset graphic" width="500" height="300" />| 
</p>
<br/><br/>

How can we tell when someone has fallen asleep? Understanding the way we fall asleep is an important problem in sleep medicine, since sleep disorders can disrupt the process of falling asleep. In the case of insomnia, subjects may fall asleep too slowly, whereas during sleep deprivation or narcolepsy, subjects fall asleep too quickly. Current methods for tracking the wake/sleep transition are time-consuming, subjective, and simplify the sleep onset process in a way that severely limits the accuracy, power, and scope of any resulting clinical metrics. We have developed a new physiologically principled method that dynamically combines information from brainwaves, muscle activity, and a novel minimally-disruptive behavioral task, to automatically create a continuous dynamic characterization of a person’s state of wakefulness.

Specifically, we propose a fully Bayesian dynamic state-space model framework for the characterization of simultaneously observed behavioral and physiological dynamics during the SOP, which is implemented using a partical filter approach. In doing so, we create a robust quantitative representation of SOP dynamics that can be used to more accurately and more precisely track the gradual transition from wakefulness to sleep.

## Video Tutorial
A video tutorial providing an overview of the general concepts underlying the code [can be found here.](https://www.youtube.com/watch?v=wAGD3Qq6n_w)
