# Credible-mathematical-model-of-the-cardiovascular-system-to-fluid-perturbation
This software code provides a mathematical model of the cardiovascular system (CVS) response to fluid perturbation, developed and validated using data collected from animal (sheep and swine) subjects. The model is built in MATLAB Simulink environment and the goal is to generate a cohort of simulated subjects against a test subject. The simulated subjects lead to a prediction envelope, which is used to evaluate predictive capability performance of the mathematical model and show its ability to generate valid virtual CVS physiological responses to fluid perturbation for the future non-clinical simulated testing setups. These testing setups may be useful as part of assessing physiological closed-loop control algorithms for automated fluid resuscitation systems. The CVS mathematical model is a low-order lumped parameter model, designed for use with virtual cohort generation tools. It takes rates of hemorrhage, urine, and fluid infusion as inputs and produces outputs for hematocrit (HCT), blood volume (BV), heart rate (HR), stroke volume (SV), cardiac output (CO), and mean arterial blood pressure (BP). The mathematical model calibration is defined based on the maximum likelihood estimation of its parameters. A compartment-based virtual cohort generation tool described in [1] is used to simulate virtual subjects and generate a prediction envelope used for model predictive capability performance. For more information about the animal study please refer to [2].

Disclaimer: This software and documentation (the "Software") were developed at the Food and Drug Administration (FDA) by employees of the Federal Government in the course of their official duties. Pursuant to Title 17, Section 105 of the United States Code, this work is not subject to copyright protection and is in the public domain. Permission is hereby granted, free of charge, to any person obtaining a copy of the Software, to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, or sell copies of the Software or derivatives, and to permit persons to whom the Software is furnished to do so. FDA assumes no responsibility whatsoever for use by other parties of the Software, its source code, documentation or compiled executables, and makes no guarantees, expressed or implied, about its quality, reliability, or any other characteristic. Further, use of this code in no way implies endorsement by the FDA or confers any advantage in regulatory decisions. Although this software can be redistributed and/or modified freely, we ask that any derivative works bear some notice that they are derived from it, and any modified versions bear some notice that they have been modified.

References

[1] Varun Kanal, Pras Pathmanathan, Jin-Oh Hahn, George Kramer, Christopher Scully, and Ramin Bighamian, Development and Validation of a Mathematical Model of Heart Rate Response to Fluid Perturbation, Scientific Reports 12, 21463 (2022).

[2] Abraham Rafie, Paul Rath, Michael Michell, Robert Kirschner, Donald Deyo, Donald Prough, James Grady, George Kramer, Hypotensive Resuscitation of Multiple Hemorrhages using Crystalloid and Colloids, Shock 22, 262-269 (2004).
