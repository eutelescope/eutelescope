This is the EUTelescope package, a generic pixel telescope data
analysis framework.

You can find all documentation online: http://eutelescope.web.cern.ch
and the source code at https://github.com/eutelescope/eutelescope

To install the development version of EUTelescope, please follow the instructions at http://eutelescope.web.cern.ch/content/installation *but* make sure that you have the git client installed on your system and download the ```eutel-git``` branch of ilcsoft-install instead of the version stated on the installation pages:
```
svn co https://svnsrv.desy.de/public/ilctools/ilcinstall/branches/eutel-git ilcinstall_eutel-git
```

Should you have any questions regarding EUTelescope, please check our
online support pages, especially the forums and the FAQ or contact us
via mail: eutelescope-coordinators@desy.de

For feature requests or bug reports please open a new issue on the issue tracker on GitHub.

Maybe you would also like to get involved in the development - all
contributions are highly welcomed!

Have you fixed the bug already or want to contribute your own processor into the main repository? Then please fork the 
project and issue a pull request using these instructions:


* create a user account on github, log in
* 'fork' the (main) project on github (using the button on the page of the main repo)
* now add the newly forked as a git remote to your EUTelescope installation and rename the original repository to 'upstream':
```
cd $ILCSOFT/Eutelescope/trunk
git remote rename origin upstream
git remote add origin https://github.com/[YOUR GITHUB USER HERE]/eutelescope
git remote -v show
```
* now edit away on your local clone! But keep in sync with the development in the upstream repository by running
```
git pull upstream master
```
on a regular basis. Replace ```master``` by the appropriate branch if you work on a separate one.
Don't forget that you can refer to issues in the main repository anytime by using the string ```eutelescope/eutelescope#XX``` in your commit messages, where 'XX' stands for the issue number, e.g.
```
[...]. this addresses issue eutelescope/eutelescope#1
```

* push the edits to origin (your fork)
```
git push origin
```
(defaults to 'git push origin master' where origin is the repo and master the branch)

* verify that your changes made it to your github fork and then click there on the 'compare & pull request' button

* summarize your changes and click on 'send' 
