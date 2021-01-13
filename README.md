Imploded cauliflower - is a Julia set for c=1/4+epsilon with epsilon >0. 

How many different ways are there to show such set ? ( program d.c creates 39 pgm images )

Here c = 0.35 



# Images
See:
* [images in commons category](https://commons.wikimedia.org/wiki/Category:Imploded_cauliflower)


## dynamic plane = z plane

![](./images/de.png "boundary using DEM/J") 
![](./images/pot.png "potential and boundary by DEM") 
![](./images/bd.png "BD/J") 
![](./images/bdb.png "boundaries of BD/J") 
![](./images/mbd.png "MBD/J") 
![](./images/mdbb.png "boundaries of MBD/J") 


![](./images/ls.png "LSM/J") 
![](./images/lc.png "boundaries of LSM/J") 
![](./images/lsc.png "LSM + boundaries of LSM/J")
![](./images/lcmbd.png "boundaries of LSM/J and MBD") 
![](./images/u.png "Unknown : boundary and slow dynamics") 
![](./images/sac.png "SAC/J + DEM/J") 
![](./images/dld.png "DLD/J + boundary by DEM") 

### procedural texture


![](./images/t0n.png "") 
![](./images/t1n.png "") 
![](./images/t2n.png "")
![](./images/t3n.png "") 
![](./images/t4n.png "") 
![](./images/t5n.png "") 
![](./images/t6n.png "") 

with boundaries 

![](./images/t5nb.png "") 
![](./images/t6nb.png "") 
![](./images/t7nb.png "") 



### Normal shading
![](./images/np.png "NP + DEM") 
![](./images/nd.png "ND + DEM") 


### blended gradients

[Description](https://github.com/adammaj1/Mandelbrot-set-with-blended-gradients)  


![](./images/blend.png "blend of normal shading and smooth escape time") 

## Inverted plane  = w plane  = 1/z plane 

![](./images/dei.png "boundary using DEM/J inv") 
![](./images/bdi.png "BD/J inverted") 
![](./images/bdbi.png "boundaries of BD/J inv") 
![](./images/sacdei.png "SAC/J + DEM/J inverted") 
![](./images/lsi.png "LSM/J inv") 
![](./images/lci.png "boundaries of LSM/J inv") 
![](./images/lsci.png "LSM + boundaries of LSM/J inv") 
![](./images/pot_i.png "potential +  boundary by DEM inverted") 

### Normal shading
![](./images/npi.png "NP + DEM inverted") 
![](./images/ndi.png "ND + DEM inverted") 
### blended gradients
![](./images/blend_i.png " inverted blend of normal shading  and smooth escape time") 


## test images

![](./images/defq.png "boundary using DEM/J and first quadrant") 
![](./images/wonz.png "W Window On Z Window") 
![](./images/zonw.png "Z Window On W Window") 




# Files
* [d.c ](./src/d.c) - c console program for creating pgm images
* [i.c](./src/i.c) - c console program for testing functions from the main program (d.c file). Prints info about point z  
* [g.sh](./src/g.sh) - bash script for conversion ( from pgm to png ) and resizing ( downscalling) using Image Magic
* [Makefile](./src/Makefile)
* [a.txt](./src/a.txt) - text output of d.c program

to run:
```
make
```

# technical notes




## Contributors

are wellcome 

See als [FF](https://fractalforums.org/programming/11/how-many-different-ways-are-there-to-show-such-set/3874) 


  
## License

[GPL](https://www.gnu.org/licenses/gpl-3.0.html)



## Git

GitLab uses:
* the Redcarpet Ruby library for [Markdown processing](https://gitlab.com/gitlab-org/gitlab-foss/blob/master/doc/user/markdown.md)
* KaTeX to render [math written with the LaTeX syntax](https://gitlab.com/gitlab-org/gitlab-foss/blob/master/doc/user/markdown.md), but [only subset](https://khan.github.io/KaTeX/function-support.html)






### Subdirectory

```git
mkdir images
git add *.png
git mv  *.png ./images
git commit -m "move"
git push -u origin master
```
then link the images:

```txt
![](./images/n.png "description") 

```

```git
git mv -f 
```




### repo



```git
cd existing_folder
git init
git remote add origin git@gitlab.com:adammajewski/implodedcauliflower.git
git add .
git commit -m "Initial commit"
git push -u origin master
```


To clone repo

```git
git clone git@gitlab.com:adammajewski/implodedcauliflower.git
```



Local repo : ~/Dokumenty/branched_ray/c/

