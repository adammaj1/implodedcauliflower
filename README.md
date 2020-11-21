Imploded cauliflower - is a Julia set for c=1/4+epsilon with epsilon >0. 

How many different ways are there to show such set ?


See:
* [images in commons category](https://commons.wikimedia.org/wiki/Category:Imploded_cauliflower)



![](./images/1.png "description") 
![](./images/2.png "description") 
![](./images/3.png "description") 
![](./images/4.png "description") 
![](./images/5.png "description") 
![](./images/6.png "description") 
![](./images/7.png "description") 
![](./images/8.png "description") 
![](./images/9.png "description") 
![](./images/10.png "description") 

![](./images/11.png "description") 
![](./images/12.png "description") 
![](./images/13.png "description") 
![](./images/14.png "description") 
![](./images/15.png "description") 
![](./images/16.png "description") 
![](./images/17.png "description") 
![](./images/18.png "description") 
![](./images/19.png "description") 
![](./images/20.png "description") 


![](./images/21.png "description") 
![](./images/22.png "description") 




# Files
* [d.c ](./src/d.c) - c console program
* [g.sh](./src/g.sh) - bash script for conversion ( from pgm to png ) and resizing ( downscalling) using Image Magic


# technical notes




## Contributors

are wellcome 


  
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

