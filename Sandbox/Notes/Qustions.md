# Questions:

1. When fitting a model, do we generate the model for the whole dataset then compare the fit of each ID to that model, or do we create a model for each ID and compare those parameters?

    - **Answer**: Curve by curve, NOT the whole dataset

2. **Starting values**:
- For the linear models, do we need to calculate starting values? 

    - **Answer**: LM for linear models, NLSLM for non-linear models.

- Pablo's code from yesterday seemed to suggest that the minimize function in python found the best starting values to minimise the residuals. If this is right then why will changing the starting values for parameter estimates make any difference to the overall fit (if minimse is finding the best parameters)?

    - **Answer**: No recommend `multistart` package in R

3. Does it matter if we use lm or nlls for linear models? And for linear models, if the data has been logged, do we need to somehow log the model if we are using the lm method with poly?

    - **Answer**: Doesnt matter
    - **Answer**: `Generate a column for logarimised data`

4. I tried to increase the maximum number of iterations when using `nlsLM()` in R, but it stopped before I could get to the maxiter specified because of some 'maxfev' threshold. I'm not really sure what this means...should/can I touch it?

    - **Answer**: maxfev: Max Func Evaluations
        - starting values
        - model
        - the data is very poor
            - `Not to worry too much!`

5. Question about the $1 \div k \times T$ scale - how to find starting values for Sharpe-Schoolfield equation - what happens when you take a log.

    - $E_{l} ~= E$
    - $E_H >> E$
