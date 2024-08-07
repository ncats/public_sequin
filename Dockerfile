FROM ncats/sctl-rshiny-complex:base-4.0.4

RUN R -e "install.packages(c('clustree', 'showtext', 'gfonts', 'patchwork'), repos = 'http://cran.r-project.org/')"

ADD install.R /tmp/
RUN R -f /tmp/install.R

COPY app.R /srv/shiny-server/
COPY global.R /srv/shiny-server/
COPY iPSCeq-server.R /srv/shiny-server/
COPY iPSCeq-ui.R /srv/shiny-server/
COPY iPSCeq-functions.R /srv/shiny-server/
COPY iPSCeq-tabs.R /srv/shiny-server/
COPY iPSCeq-pack-load.R /srv/shiny-server/
COPY data/ /srv/shiny-server/data/
COPY markdown/ /srv/shiny-server/markdown/
COPY www/ /srv/shiny-server/www
COPY iPSCeq-modules.R /srv/shiny-server/

COPY shiny-server.conf /etc/shiny-server/shiny-server.conf
