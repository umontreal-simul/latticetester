# Updating LatticeTester documentation	

 To update the documentation of *LatticeTester* , run the following commands on your clone of the LatticeTester repository:	

 ```	
./waf configure --build-docs	
./waf build	
rm -r doc/html	
git clone https://github.com/umontreal-simul/latticetester.git doc/html	
cd doc/html	
git checkout -b gh-pages origin/gh-pages	
rm -r !('README.md')
cp ../../build/doc/html/* .
git status	
git add .	
git commit -m "update documentation"	
git push	
```	

 These commands should compile the source code and the documentation, commit and push the documentation changes to the branch `gh-pages` of *LatticeTester*.
