FROM    nfcore/base:1.9
LABEL authors = 'Francesc Català-Moll' \
      dasdfsa = 'Docker image containing all software requeriments for the nfcovid pipeline' \
      date = '2020-09-06'

# Install conda envirorment 
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a 

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nfcovid-1.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nfcovid-1.0 > nfcovid-1.0.yml

# Instruct R processes to use these empty files instead of clashing with a local version
RUN touch .Rprofile
RUN touch .Renviron

