# Water Molecule Counting in Lipid Membranes

This script counts the number of water molecules between the upper and lower halves of a lipid membrane in each frame of a molecular dynamics (MD) simulation trajectory.

## Description

The script uses the [MDAnalysis](https://www.mdanalysis.org/) to process the trajectory and calculate the number of water molecules located between the upper and lower halves of the membrane for each frame in the trajectory. 
It also provides options for saving residue IDs, calculating mean and error, and plotting the data.

Originally, this script is used for calculating water molecules inside the hydrophobic core of a lipid membrane. To estimate the number of water molecules in the hydrophobic core (in the region of lipid tails), 
we used the following simple method. We defined the boundaries of the hydrophobic core using the C32 and C22 carbon atoms of lipids (first C-atoms after carboxyl group for each lipid). 
For each frame of the MD trajectory, the minimum distance r between C-atoms of different half-layers is determined. Then we select molecules which distances to each half-layer is less than r. 
Indeed, this is only possible for those water molecules that are located within these boundaries. 
All water molecules outside the membrane will have a greater distance than r to one or another layer of C-atoms.  

This method was first described and applied in [this article](https://doi.org/10.1016/j.molliq.2024.123948). Please cite:  
>Yakush, E.A., Shelepova, E.A. and Medvedev, N.N., 2024. Mechanism of water transport through the lipid membrane with trichogin GA IV. Molecular dynamics study. Journal of Molecular Liquids, 396, p.123948.


## Installation

1. Clone the repository:

   ```sh
   git clone https://github.com/ElAsh1906/water-count.git
   cd water-count
   ```

2. Install the required Python packages:  
   
   ```sh
   pip install MDAnalysis pandas numpy matplotlib seaborn tqdm
   ```
   
3. For installation as a utility, you can use the following (Note: You may need to provide your password for sudo.):

   ```sh
   chmod +x wcount.py
   sudo mv wcount.py /usr/local/bin/wcount
   source ~/.bashrc
   ```

## Usage
To run the script, use the following command:

`
wcount [-h] [-c <.gro>] [-s <.tpr>] [-f <.xtc>] [-o <.dat/.txt/.csv/...>] [-l <lipid>] [-la <atom name> [<atom name> ...]] [-w <water>] [-wa <atom>] [-d <distance>] [-rw <rolling window>] 
[-bf <first frame>] [-b <start time>] [-ef <last frame>] [-e <last time>] [-step <step>] [-ts <time step>] [--save_resids] [--calc] [--plot]
`

#### Options to specify input files are:

`-c` &emsp;&emsp; <.gro> &emsp;&emsp; (conf.gro)  
&emsp;&emsp;&emsp;&emsp;Input topology file.  
`-s` &emsp;&emsp; <.tpr> &emsp;&emsp; (topol.tpr)  
 &emsp;&emsp;&emsp;&emsp;Input structure file (optional).  
`-f` &emsp;&emsp; <.xtc> &emsp;&emsp; (traj_comp.xtc)  
&emsp;&emsp;&emsp;&emsp;Input trajectory file.  

#### Options to specify output files are:

`-o` &emsp;&emsp; <.dat/.txt/.csv/...> &emsp;&emsp; (water_mols.csv)  
&emsp;&emsp;&emsp;&emsp;Output file for results.  
`--save_resids` &emsp;&emsp; (resids.csv)  
&emsp;&emsp;&emsp;&emsp;Flag to save all residue IDs. 

#### Other options are:

`-l, --lipids` &emsp;&emsp; <lipid> &emsp;&emsp; (DOPC)  
&emsp;&emsp;&emsp;&emsp;Lipid model used in the simulation.  
`-la` &emsp;&emsp; <atom name>  
&emsp;&emsp;&emsp;&emsp;Atoms in the lipid molecule to use for analysis.  
`-w` &emsp;&emsp; <water> &emsp;&emsp; (SOL)  
&emsp;&emsp;&emsp;&emsp;Residue name of the molecule to be counted.  
`-wa` &emsp;&emsp; <atom> &emsp;&emsp; (OW)  
&emsp;&emsp;&emsp;&emsp;Main atom name of the molecule to be counted.  
`-d` &emsp;&emsp; <distance> &emsp;&emsp; (0 Å)  
&emsp;&emsp;&emsp;&emsp;Distance threshold for counting molecules in Å.  
`-rw` &emsp;&emsp; <rolling window> &emsp;&emsp; (10)  
&emsp;&emsp;&emsp;&emsp;Size of the rolling window for smoothing the data.  
`-bf` &emsp;&emsp; <first frame> &emsp;&emsp; (0)  
&emsp;&emsp;&emsp;&emsp;Frame number to start analysis from.  
`-b` &emsp;&emsp; <start time> &emsp;&emsp; (0 ps)  
&emsp;&emsp;&emsp;&emsp;Start time for analysis in ps.  
`-ef` &emsp;&emsp; <last frame> &emsp;&emsp; <-1>  
&emsp;&emsp;&emsp;&emsp;Frame number to end analysis at.  
`-e` &emsp;&emsp; <last time> &emsp;&emsp; <-1>    
&emsp;&emsp;&emsp;&emsp;End time for analysis in ps.  
`-step` &emsp;&emsp; <step> &emsp;&emsp; (1)  
&emsp;&emsp;&emsp;&emsp;Interval for frame-based analysis.  
`-ts` &emsp;&emsp; <time step>  
&emsp;&emsp;&emsp;&emsp;Time step interval for analysis.   
`--calc`  
&emsp;&emsp;&emsp;&emsp;Flag to calculate mean and error of the results.  
`--plot`  
&emsp;&emsp;&emsp;&emsp;Flag to plot the data.  

## Output

The output file water_mols.csv will contain the number of water molecules between the upper and lower halves of the membrane for each frame.
If the --save_resids flag is used, a file resids.csv will be created with the residue IDs of the water molecules for each frame. 

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Contact

If you have any questions or need further assistance, feel free to contact elenayakush43@gmail.com.
















