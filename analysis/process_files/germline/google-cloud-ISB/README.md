Preliminary steps (any order):

* Build dockerfile and push to repository of choice

* Install dsub

  see [https://github.com/googlegenomics/dsub](https://github.com/googlegenomics/dsub)

* Make lists of samples ids to process (if needed)

   ./make_lists.sh

    This action creates the directory analysisID_lists/ and files therein



Main steps

1. In script annotate.sh, check for the appropriate paths and filenames, dsub location, docker image repository location, as well as cloud project information and container parameters, etc.

2. Activate dsub.

3. Run

   ./annotate.sh  $cancerType

   This action typically launches several hundred jobs, which is the recommended approach.
