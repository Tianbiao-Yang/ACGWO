# ACGWO
An improved chaotic gray wolf optimization algorithm was used
1) Improvement and Application of GWO Algorithm
Grey wolves prey on the smaller animals in nature, relying on their leadership hierarchy and hunting mechanism. In 2014, GWO algorithm is applied to achieve the simulation by building a 4-layer pyramid hierarchical management system, which performs the social hunting behavior. The management system is divided into α (the first layer), β (the second layer), δ (the third layer), ω (the fourth layer). During the hunt, α, β and δ perform hunting and search behaviors, ω are responsible for the safety and integrity in the dominance structure, and ultimately complete the predation task. Further, α is guider, who mostly responsible for making decisions. Nevertheless, β are subordinate wolves that help α in collective activities, which is the best candidate to be the leader. δ take orders from α and β, but they dominate the ω.
Grey wolves recognize and encircle prey, and the location of wolf (X, Y) is usually updated according to the position of prey (X', Y'). In the population, the best individual is α, which defined as the optimal solution in the history. The fitness values ranking second (suboptimal solution) and third (final solution) are denoted as β and δ, while the remaining individuals are denoted as ω. In the d-dimensional search space, the number of gray wolves is assumed to be N, the number of iterations is expressed by t, and the position vector of the i-th gray wolf is calculated by the formula (S1, S2). 
                                  (S1)
                               (S2)
Where   represents the position vector of the i-th gray wolf,   represents the prey’ position vector,   is the swing factor,   is the direction,   is the convergence factor,   and   represent the random vectors in [0,1]. Therefore, the position vector of the prey can be calculated as formula (S3).
          (S3)
After each iteration by GWO algorithm, the three best solutions are assigned to α, β, δ. Then, the position of other wolves (search agents) is updated according to the position of α, β and δ. In other words, α, β and δ estimate the position of the prey, and the other wolves randomly update their position around the prey. Finally, the position of the prey is the optimal solution after multiple iterations.
2) The Population Initialization Based on CT
To maintain the diversity of the population and make the initial population uniformly distributed, the chaotic algorithm was introduced in GWO algorithm to accelerate the convergence speed of the algorithm. The chaos variables are firstly linearly mapped into the region of the optimization variables, and chaos transform is then utilized to search. Chaos randomness and ergodicity can avoid the tendencies of falling into local minimal values during the search process, thus overcoming the shortcomings of GWO algorithm. By studying the chaotic map sequence and the traditional optimization algorithm, the Logistic map with good ergodicity and uniformity is selected. The mapping equation is described as formula (S4).
             (S4)
The population is in a chaotic state when µ = 4, and its output is equivalent to a random variable between [0,1], which can traverse each number of interval [0,1] without any repetition. The mapping structure is simple with good traversal uniformity and better iterative speed, making the objective function easier to converge. As observed in Fig.S2, the Chaotic Grey Wolf Optimization (CGWO) algorithm could be clearly explained. 
 
Fig. S II1 CGWO algorithm

3) Improved Convergence Factor 
The coordinated problem between global search and local search were common in the group intelligence algorithms. The strong global search ability can ensure the diversity of the population, and the strong local search capability can guarantee the precision of results. Therefore, it is important to deal with the balance between the global search capability and the local search capability in the GWO algorithm. The convergence factor a in the basic GWO algorithm is decreases linearly as the number of iterations increases, but the algorithm does not converge linearly during the process of continuous convergence. The linear convergence of the convergence factor a can’t fully reflect the actual optimization process [32]. To maintain the balance between global search and local search, a nonlinear convergence factor was introduced:
                       (S5)
Where: e represents the number of natural logarithms, t represents the actual number of the current iterations, m represents the maximum number of iterations.  
As seen in Fg. S3, the convergence factor a is nonlinearly decreases with the number of iterations from 2 to 0. At first, the attenuation of the initial convergence factor is reduced so as to find the global optimal solution as soon as possible. Lately the attenuation of the convergence factor is increased so that the local optimal solution can be accurately searched. Therefore, the nonlinear convergence method can maintain the balance between the diversity of global search and the accuracy of local search so as to improve the stability of the convergence. Here, the nonlinear convergence factor a could also be introduced in the basic GWO algorithm (AGWO). 
 
Fig. S II2  The nonlinear convergence factor a is introduced in AGWO algorithm
 As observed in Fig. 1, the optimization steps of ACGWO algorithm are explained. Firstly, the algorithm parameters were set. Secondly, the fitness values of individuals in the population are calculated and sorted, then the first three values with better fitness are assigned to α, β and δ. Thirdly, the position of ω and the chaotic sequence are updated. Fourthly, the number of iterations and convergence factor a were updated, to calculate the value of the improved   and update the parameters ( , , , , and xn) etc. Finally, the optimal solution (α) could be obtained. 
