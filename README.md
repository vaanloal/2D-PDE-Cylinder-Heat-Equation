# 2D-PDE-Cylinder-Heat-Equation
Utilization of finite differences - method of lines (MOL) - to solve an unsteady state 2D Partial Differential Equation of an cork cylinder, in MATLAB.

to write...

# Model

The dimensionalizes model can be written as:

![Dimensionalizes model.](/images/dimensionalized_model.png)

At the centerline of cylinde the model becomes the next equation. At r = 0, because BC and the term 1/r the model has one indetermination. Applyg L'Hôpital's rule we can arrive at:

![dimensionalized_model_at_center.png](/images/dimensionalized_model_at_center.png)

The Initial and Bounday Conditions are given by the following equations:

![dimensionalized_model_BC.png](/images/dimensionalized_model_BC.png)
