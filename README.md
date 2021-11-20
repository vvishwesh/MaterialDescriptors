# MaterialDescriptors

The software calculates a fixed length set of descriptors for a given chemical composition. The descriptors are derived from elemental compositions include, for instance, indices such as the stoichiometric average of the elemental properties (e.g. atomic weight, volume, density, electronegativity, free orbital electrons etc). These attributes capture the fraction of the elements present. Other descriptors focus on values such as the maximum, minimum, fraction-weighted mean, average deviation, and mode (i.e. the property of the most prevalent element) of the attributes.

## Compilation
To compile the program

```
 cd build
 bash build.sh
```
The above commands will create ```MaterialDescriptors.jar``` which you can use to calculate the descriptors.


## Usage
To use the program

```
Usage:java -jar MaterialDescriptors.jar inputFile header (1/0/yes/no/true/false) outputFile

inputfile: a text file containing the chemical compositions.
header: whether the input file contains a header line
outputfile: name of the file to which the descriptors will be written.
```
