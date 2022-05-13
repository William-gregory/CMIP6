# CMIP6
The code in this repository can be used to download CMIP6 model outputs, and subsequently run complex networks analysis, as outlined in the article [Gregory et al., 2022](https://tc.copernicus.org/articles/16/1653/2022/tc-16-1653-2022.pdf).

<img src="https://github.com/William-gregory/CMIP6/blob/main/images/network_connections.png" width="600" height="500">

# Downloading data
The executable file `cmip6_download_and_extract.py` must be run on the UK computer server JASMIN, where it will search the JASMIN directories for available model ensembles which match the input description provided by the user. If it cannot find any files, it will resort to downloading the files remotely from the CMIP6 data portal (https://esgf-node.llnl.gov/projects/cmip6/) by using the download script created by Thiago Loureiro (https://github.com/tloureiro/cmip6_downloader). The extracted variable(s) are then regridded to a 50km polar stereographic grid and all ensembles for that model are then saved in a pickle file.

Note that the files `cmip6_download_and_extract.py` and `cmip6_downloader.py` are assumed to be in the same directory

Here are two examples. The first shows an example of where the data are available on JASMIN:

![alt text](https://github.com/William-gregory/CMIP6/blob/main/images/example1.png)

The second shows how the data are unavailable on JASMIN and must be downloaded remotely:

![alt text](https://github.com/William-gregory/CMIP6/blob/main/images/example2.png)

It is also worth noting that sometimes the remote download can fail for various reasons, and in some cases can be overcome simply by re-running the code. For specific bugs relating to the remote download, please contact Thiago Loureiro (thiago@tloureiro.com).


