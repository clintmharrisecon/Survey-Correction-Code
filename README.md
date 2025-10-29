# Code
These files illustrate how to use data request and response timing to construct instruments for nonresponse corrections as proposed by Harris, Eckhardt, and Goldfarb, 2025. The paper can be found at https://arxiv.org/abs/2404.17693.

There are 3 files.

1. Selection_In_Surveys_Example.do
Generates a synthetic version of the Norway in Corona Times survey using estimates reported in "Selection in Surveys: Using Randomized Incentives to Detect and Account for Nonresponse Bias" available at https://www.nber.org/papers/w29549.
Then performs nonresponse corrections for both continuous and binary variables using request instruments.

2. Simulated_Example_continuous.do
Generates a selectively observed continuous survey variable and performs a nonresponse correction using request instruments, as well as diagnostic tests.

3. Simulated_Example_binary.do
Generates a selectively observed binary survey variable and performs a nonresponse correction using request instruments, as well as diagnostic tests.


