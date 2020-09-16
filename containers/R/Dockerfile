# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/pandora1_paper_r:pandora_paper_tag1

# This container has R and the libraries needed to build plots in this repo
# to build: sudo docker build . -t leandroishilima/pandora1_paper_r:pandora_paper_tag1
FROM rocker/tidyverse:4.0.0

ARG extra_packages="c(\"gridExtra\")"
RUN echo "install.packages($extra_packages, repos=\"https://cran.rstudio.com\")" | R --no-save
RUN cp /usr/bin/python3 /usr/bin/python
