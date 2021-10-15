# Shiny-server

This server offers user GUI to select .count files and generate gene expression file.

## Installation

```bash
sudo apt-get install r-base
sudo apt-get install gdebi-core
sudo gdebi shiny-server-<version>.deb
```

## Useful commands
```bash
sudo systemctl start shiny-server
sudo systemctl stop shiny-server
sudo systemctl restart shiny-server
```
## Configuration files

```bash
/etc/shiny-server/shiny-server.conf 
```
Basic port number : listen 3838;  
Basic site directory : site_dir /srv/shiny-server;  
Basic site directory : log_dir /var/log/shiny-server;  

## Manual
Shiny server.docs

## License
[MIT](https://choosealicense.com/licenses/mit/)
