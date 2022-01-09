## Cite Paper
Karve, P. M., Guo, Y., Kapusuzoglu, B., Mahadevan, S., & Haile, M. A. (2020). Digital twin approach for damage-tolerant mission planning under uncertainty. Engineering Fracture Mechanics, 225, 106766.

Please, cite this repository using: 

@article{karve2020digital,
  title={Digital twin approach for damage-tolerant mission planning under uncertainty},
  author={Karve, Pranav M and Guo, Yulin and Kapusuzoglu, Berkcan and Mahadevan, Sankaran and Haile, Mulugeta A},
  journal={Engineering Fracture Mechanics},
  volume={225},
  pages={106766},
  year={2020},
  doi={10.1016/j.engfracmech.2019.106766},
  publisher={Elsevier}
}
	
  
# Digital-twin-approach-for-damage-tolerant-mission-planning-under-uncertainty
The digital twin paradigm that integrates the information obtained from sensor data, physics models, as well as operational and inspection/maintenance/repair history of a system (or a component) of interest, can potentially be used to optimize operational parameters of the system in order to achieve a desired performance or reliability goal. In this article, we develop a methodology for intelligent mission planning using the digital twin approach, with the objective of performing the required work while meeting the damage tolerance requirement. The proposed approach has three components: damage diagnosis, damage prognosis, and mission optimization. All three components are affected by uncertainty regarding system properties, operational parameters, loading and environment, as well as uncertainties in sensor data and prediction models. Therefore the proposed methodology includes the quantification of the uncertainty in diagnosis, prognosis, and optimization, considering both aleatory and epistemic uncertainty sources. We discuss an illustrative fatigue crack growth experiment to demonstrate the methodology for a simple mechanical component, and build a digital twin for the component. Using a laboratory experiment that utilizes the digital twin, we show how the trio of probabilistic diagnosis, prognosis, and mission planning can be used in conjunction with the digital twin of the component of interest to optimize the crack growth over single or multiple missions of fatigue loading, thus optimizing the interval between successive inspection, maintenance, and repair actions.

## RBDO_Uncert_a0_C_m
The folder contains optimization of load history and number of cycles under uncertainty.
The initial crack size, and the model parameters are uncertain parameters that are quantified by using the Metropolis-Hastings algorithm which is a Markov chain Monte Carlo method.
Gaussian process surrogate model is constructed to repsent the fatigue crack growth model.
