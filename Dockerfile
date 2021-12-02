FROM continuumio/miniconda3

WORKDIR /app

# Create the environment:
COPY polyfun.yml .
RUN conda env create -f polyfun.yml
RUN chmod 775 pipeline.sh
# Make RUN commands use the new environment:
RUN echo "conda activate polyfun" >> ~/.bashrc
#SHELL ["/bin/bash", "--login", "-c"]


#RUN R -e "install.packages('susieR')"

# The code to run when container is started:

#COPY . ./
ENTRYPOINT ["./entrypoint.sh"]
