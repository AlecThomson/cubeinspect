# cubeinspect

A tool for quickly inspecting channels in large astronomical data cubes, stored in FITS. Written for use by the POSSUM collaboration.


## Installation

Use the package manager [pip](https://pypi.org/) to install cubeinspect.

```bash
pip install cubeinspect
```

## Usage

```
$ cubei -h
usage: cubei [-h] [-v] [-t] [-c vmin vmax] [-d DPI] fitsfile channel

    Open big file and make an average image of a selected channel.
    Image will be squared so that Stokes Q and U look reasonable.

    Saves a png with same name as FITS image + .medimage.png



positional arguments:
  fitsfile      FITS file to open
  channel       Channel to inspect

optional arguments:
  -h, --help    show this help message and exit
  -v            Verbosity.
  -t            Make thumbnail (takes extra time!)
  -c vmin vmax  Limits for image (defaults to vmax=std(image))
  -d DPI        Resolution of saved images in DPI.
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
[BSD 3-clause](https://choosealicense.com/licenses/bsd-3-clause/)