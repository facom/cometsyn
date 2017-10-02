```
//////////////////////////////////////////////////////////////////////////////////
//    ___                     _   __             
//   / __\___  _ __ ___   ___| |_/ _\_   _ _ __  
//  / /  / _ \| '_ ` _ \ / _ \ __\ \| | | | '_ \ 
// / /__| (_) | | | | | |  __/ |__\ \ |_| | | | |
// \____/\___/|_| |_| |_|\___|\__\__/\__, |_| |_|
//                                   |___/       
// 2013,2017 [)] Jorge Zuluaga, jorge.zuluaga@udea.edu.co
// Instituto de Fisica / FCEN - Universidad de Antioquia
//////////////////////////////////////////////////////////////////////////////////
```

Presentation
------------

Simulates a Comet.

Getting a copy
--------------

To get a copy of the newest version of this project just execute:

```
$ git clone --branch BHMcalc2 http://github.com/facom/BHMcalc.git
```

For the oldest version just remove the `--branch` option.

Instructions for the contirbutor
--------------------------------

1. Generate a public key of your account at the client where you will
   develop contributions:
   
   ```
   $ ssh-keygen -t rsa -C "user@email"
   ```

2. Upload public key to the github FACom repository (only authorized
   for the FACom repository manager), https://github.com/facom.

3. Configure git at the client:

   ```
   $ git config --global user.name "Your Name"
   $ git config --global user.email "your@email"
   ```

4. Get an authorized clone of the project:

   ```
   $ git clone git@github.com:facom/BHMcalc.git
   ```

5. Checkout the branch you are interested in (e.g. BHMcalc2):

   ```
   $ git checkout -b BHMcalc2 origin/BHMcalc2
   ```

6. Checkout back into the master:

   ```
   $ git checkout master
   ```

Package structure
-----------------

* cometsyn.cpp: main routines.

* simulate.cpp: comet simulation routine.


Running a simulation
--------------------