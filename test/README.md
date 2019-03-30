# Unit Tests

It is crucial that this fork of combine is well tested and validated such that it can be trusted for physics analysis.

Therefore, I encourage everyone who wants to use it to test it to validate against the official version of combine and implement this validation as a unit test here. This way, the next user who has the same usecase does not to do the validation herself.

Since the reference output of the original combine should not change, it is sufficient to extrast the relevant output information of a given method in combine from the output ROOT file and hardcode it in these unit tests.

The unit tests require an agreement of __at least 6 digits__. Please don't require any unit test with looser requirements!

Travis CI is used to automatically run the unit tests for the validation of every commit. The testing platform has the following environment:

* Ubuntu 16.04 Xenial Xerus
* Python 2.7
* ROOT 6.16

## List of implemented validations

You can fine [a list of the equivalent commands in the original combine](validations.txt) which got implemented as unit tests.

The reference version of combine was taken on 23. March 2019 from the `81x-root606` branch at the commit with the hash `5cc169efd9233011924e5fbd468ef05be044ed39`.
