# Evaluating stochastic seeding strategies in networks

This contains code for conducting the simulations and analyses in the paper:
Chin, A., Eckles, D., &amp; Ugander, J. Evaluating stochastic seeding strategies in networks.
> When trying to maximize the adoption of a behavior in a population connected by a social network, it is common to strategize about where in the network to seed the behavior, often with an element of randomness. Selecting seeds uniformly at random is a basic but compelling strategy in that it distributes seeds broadly throughout the network. A more sophisticated stochastic strategy, one-hop targeting, is to select random network neighbors of random individuals; this exploits a version of the friendship paradox, whereby the friend of a random individual is expected to have more friends than a random individual, with the hope that seeding a behavior at more connected individuals leads to more adoption. Many seeding strategies have been proposed, but empirical evaluations have demanded large field experiments designed specifically for this purpose and have yielded relatively imprecise comparisons of strategies. Here we show how stochastic seeding strategies can be evaluated more efficiently in such experiments, how they can be evaluated "off-policy" using existing data arising from experiments designed for other purposes, and how to design more efficient experiments. In particular, we consider contrasts between stochastic seeding strategies and analyze nonparametric estimators adapted from policy evaluation and importance sampling. We use simulations on real networks to show that the proposed estimators and designs can increase precision while yielding valid inference. We then apply our proposed estimators to two field experiments, one that assigned households to an intensive marketing intervention and one that assigned students to an anti-bullying intervention.  

Some of this depends on data sets not included here:
- The data from Paluck, Shepherd &amp; Aronow (2016) about the anti-bullying intervention and the associated social networks is not included. A version of this data can be acquired from ICPSR [here](https://www.icpsr.umich.edu/web/civicleads/studies/37070/). However, our analysis was conducted using a different version of this data we received from the authors under a different data use agreement. From inspection of the ICPSR version, we believe our analyses should be possible to repeat using that data, but additional data preparation will be required to use the provided code with the version of this data on ICPSR.
- The data from Cai, De Janvry &amp; Sadoulet (2015) about the marketing intervention in Chinese villages is not included. It is readily available as part of that paper's replication materials [here](http://doi.org/10.3886/E113593V1). 
