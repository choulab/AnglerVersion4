# AnglerVersion4
AnglerVersion4 enables the design of HCR 3.0 smFISH probes targeted to single RNAs with homologous counterparts.

## Previous versions
The previous version of this code (v.1.0.0) can be downloaded [here](https://github.com/choulab/AnglerVersion4/releases).

## Usage
The current version is best used with [Docker](https://www.docker.com/), which will take care of installing BLAST, Muscle, and NUPACK and all python dependencies for you in an isolated environment. Once Docker is installed on your system, you can run `docker compose build` from the repository root directoy to build the angler image. Then you can run angler in a Docker container with `docker run`, mounting the input and output directories into the container as necessary and passing along any of the other arguments that the script accepts, which can be listed by running `docker run --rm angler4 --help`.

There is also a convenience shell script that wraps the docker command and simplifies the process. It's called [`run-angler.sh`](./run-angler.sh). Here is an example usage:

```bash
EXAMPLE: ./run-angler.sh \
  --input_dir /home/my-user/AnglerVersion4/app/static/mrna_fasta \
  --output_dir /home/my-user/angler-results-tmp \
  --log_dir /home/my-user/AnglerVersion4/app/log \
  --debug
```

The above will run the command in a docker container and mount the input and output directories so the files and results are available on your machine (and not just inside the container). **Note**: take care to ensure that both directories exist on your machine, otherwise Docker will create them and, with the input directory empty, the command will fail!

If you don't wish to use Docker, you will need to install all the necessary dependencies on your machine. For a list of these dependencies see the [Dockerfile](./Dockerfile), which includes scripts for installing them in a unix-like environment and can be copied or adapted for your operating system.
