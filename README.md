
2023.11.19

CCEX ![](https://komarev.com/ghpvc/?username=your-github-username&style=flat-square) ![](https://img.shields.io/github/watchers/{username}/{repo-name}.svg)
===============

![License](https://img.shields.io/badge/License-CQML-blue?style=flat-square&logo=c&logoColor=white&labelColor=e28a2b)<br/>
*This code is for the simulation of qubit coherence dynamics<br/>

**Document** [How to use "CCEX"](https://cce-x-refactoring.vercel.app/index.html)

**TODO**
* Add additional spins for spin defects
* Bath tensor file    
* qintmap self, easy for zfs tensor -> { D= xxx , E= xxx }
* pulse operator for each qubit
<br/>

**Useful page**
* [How to write "README.md"](https://docs.github.com/ko/get-started/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax)

**Used language**<br/>
![C](https://img.shields.io/badge/c-%2300599C.svg?style=for-the-badge&logo=c&logoColor=white)
![C++](https://img.shields.io/badge/c++-%2300599C.svg?style=for-the-badge&logo=c%2B%2B&logoColor=white)
![GitHub](https://img.shields.io/badge/github-%23121011.svg?style=for-the-badge&logo=github&logoColor=white)
![Git](https://img.shields.io/badge/git-%23F05033.svg?style=for-the-badge&logo=git&logoColor=white)
<br/>

**External library**<br/>
Eigen (For details, see [here](https://eigen.tuxfamily.org/index.php?title=Main_Page))<br/>
uthash (For details, see [here](https://troydhanson.github.io/uthash/userguide.html#_a_hash_in_c))<br/>
<br/>

> [!NOTE]
> Highlights information that users should take into account, even when skimming.

> [!IMPORTANT]
> intel MPI library is required

> [!WARNING]
> Critical content demanding immediate user attention due to potential risks.



# How to install and connect the library of mpich (mpi stuff) in macos
#### mpich
brew install mpich
mpicxx --version

$ mpicxx -v
  '-I /opt/homebrew/Cellar/mpich/4.2.1/include' [-Wunused-command-line-argument]
  The include directory : /opt/homebrew/Cellar/mpich/4.2.1/include'

