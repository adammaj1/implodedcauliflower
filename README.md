Imploded cauliflower - is a Julia set for c=1/4+epsilon with epsilon >0. 

How many different ways are there to show such set ?

Here c = 0.35 



# Images
See:
* [images in commons category](https://commons.wikimedia.org/wiki/Category:Imploded_cauliflower)


## dynamic plane = z plane

![](./images/1.png "boundary using DEM/J") 
![](./images/2.png "BD/J") 
![](./images/3.png "boundaries of BD/J") 
![](./images/4.png "MBD/J") 
![](./images/5.png "boundaries of MBD/J") 


![](./images/6.png "LSM/J") 
![](./images/7.png "boundaries of LSM/J") 
![](./images/8.png "LSM + boundaries of LSM/J")
![](./images/9.png "boundaries of LSM/J and MBD") 
![](./images/10.png "Unknown : boundary and slow dynamics") 
![](./images/11.png "SAC/J + DEM/J") 
![](./images/12.png "DLD/J + boundary by DEM") 


## Inverted plane  = w plane  = 1/z plane 

![](./images/13.png "boundary using DEM/J inv") 
![](./images/14.png "BD/J inverted") 
![](./images/15.png "boundaries of BD/J inv") 
![](./images/16.png "SAC/J + DEM/J inverted") 
![](./images/17.png "LSM/J inv") 
![](./images/18.png "boundaries of LSM/J inv") 
![](./images/19.png "LSM + boundaries of LSM/J inv") 

## test images

![](./images/20.png "boundary using DEM/J and first quadrant") 
![](./images/21.png "W Window On Z Window") 
![](./images/22.png "Z Window On W Window") 




# Files
* [d.c ](./src/d.c) - c console program
* [g.sh](./src/g.sh) - bash script for conversion ( from pgm to png ) and resizing ( downscalling) using Image Magic


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
gitm mv -f 
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

