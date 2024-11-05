BESTDR is an R package that uses Stan to estimate parameters for continuous time markov
branching processes approximated from a normal distribution. The package requires a user
to specify the model based on the transition probabilities of types along with a possible
"dose-response" statistical model connecting each rate to an underlying set of parameters.
The model is translated into Stan code that can be compiled and run.