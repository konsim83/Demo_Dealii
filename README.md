# C++ Example for Deal.ii


This C++ code is used for teaching purposes and is a modification of
the
[step-7](https://www.dealii.org/current/doxygen/deal.II/step_7.html)
example program of the Deal.ii library.


*You will need:*

 * A **Linux** distribution (I recommend Ubuntu 20.04 for beginners)
 * A **C/C++ compiler** (see package repositories)
 * **cmake** v2.8.12 or higher
 * A working installation of **[Deal.ii](www.dealii.org)** v9.1.1 or higher
 * **[Paraview](www.paraview.org)** for the visualization (free of charge)

---
**HINT**
The Deal.ii library v9.1.1 is available in the package repositories of
Ubuntu 20.04. This way you do not need to compile it yourself.	
In case you need to compile it follow the instructions linked
[here](https://www.dealii.org/current/index.html). You do not need to
install all dependencies -- a minimum is sufficient here.	
The manual of the library is available
[here](https://www.dealii.org/9.1.1/doxygen/deal.II/index.html). A
good number of video tutorials can be found
[here](https://www.math.colostate.edu/~bangerth/videos.html). They are
for an earlier release so the syntax might vary a bit but the
principles stayed essentially the same.

---

*Optional:*

* **[Eclipse IDE](https://www.eclipse.org/)** (free of charge)
  or any IDE that you prefer. Note that cmake may be able to
  generate project files also for your IDE

---


### Building and running the program

To build the project together with Eclipse project files you must first clone the repository:

```
git clone https://github.com/konsim83/Demo_Dealii.git
```
We want an out-of-source-build with build files in a folder parallel to the code:

```
mkdir build-Demo_Dealii
cd build-Demo_Dealii
```
Then create the build files with `cmake` (including project files for Eclipse):

```
cmake -DDEAL_II_DIR=/path/to/dealii -DCMAKE_ECLIPSE_MAKE_ARGUMENTS=-jN -G"Eclipse CDT4 - Unix Makefiles" ../Demo_Dealii
```
where `N` is the number of cores on your machine. You can now import an existing project in Eclipse (Eclipse).

To generate the executable in debug mode type

```
make debug
make -jN
```
If you want to produce a faster reslease version type

```
make release
make -jN
```
To run the executable with an example parameter file run

```
./helmholtz
```
Open the `vtk` files with Paraview to see the result. 
