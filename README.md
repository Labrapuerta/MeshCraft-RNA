# MeshCraft-RNA
## A Specialized Pipeline for RNA-Amino Acid Binding Meshes
Meshcraft RNA is a cutting-edge computational pipeline designed to generate meshes of amino acid residues that are capable of binding to RNA. This pipeline is an integral component of our methodology, serving as the bridge between computational modeling and practical application in RNA-aptamer interactions. 

It integrates GMLSite, created by Yidong Song and Yuedong Yang (2003), a tool designed for the identification of ligand-binding sites in protein structures. GLMsite uses Generalized Linear Models (GLMs) to predict the likelihood of interaction between amino acids and ligands. The amino acid residues identified by GLMsite are directly fed into MeshCraft RNA for mesh generation and electric potential mapping.

Meshcraft RNA also incorporates Blender, an open-source tool used for creating, editing, and transforming 3D models. With Blender, the number of faces on each mesh are reduced, making them computationally more manageable, and ensuring that the models are watertight.

## Currently we are dealing with compatibility problems of Docker and Nucleic acid binding, to solve this temporarily we developed this google colab to mesh the processed output of the program

### Get the PDB file and it's aminoacid score using 

- [Nucleic acid binding](https://github.com/biomed-AI/nucleic-acid-binding/tree/main) 

### Upload the files to the google colab notebook and follow the instructions

- [Colab Notebook](https://colab.research.google.com/drive/1EYucZ4VVxXrWLzk4PZJi365-zcfUlQsO?usp=sharing)

Our colab notebook lets you visualize the results and scores for each group
![image](https://github.com/Labrapuerta/MeshCraft-RNA/blob/main/colab%20Mesh-craft.png)

## Prerequisites
The program is containerized using Docker for easy execution.

### Install Docker Desktop
Before you can run the program, you'll need to have Docker Desktop installed on your system.
- [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/)
- [Docker Desktop for Linux](https://docs.docker.com/desktop/install/linux-install/)
- [Docker Desktop for macOS](https://docs.docker.com/desktop/install/mac-install/)

## Getting Started
1. **Clone the Repository**:

   ```bash
   git clone https://github.com/Labrapuerta/MeshCraft-RNA.git
   ```

2. **Build the Docker Image**:

   Navigate to the project directory and run the following command to build the Docker image. This image contains all the necessary dependencies for the program.

   ```bash
   docker build . -t mashcraft-rna
   ```
3. **Run the Docker Container**:

   After building the image, you can start the Docker container. Make sure Docker Desktop is running on your system. Execute the following command:

   ```bash
   docker run -p 8080:8080 mashcraft-rna
   ```
4. **Access Code Server**:

   Once the container is running, you can access Code Server by opening a web browser and navigating to [http://localhost:8080/](http://localhost:8080/). 
    You will be asked to enter a password, to find the system-generated password, follow these steps:
   
    1. Open a new terminal window.
    2. Run the following command to access the Docker container:

        ```bash
        docker exec -it mashcraft-rna /bin/bash
        ```

    3. Once inside the container, use the following command to view the configuration file:

         ```bash
         cat /root/.config/code-server/config.yaml
         ```

    4. In the output, look for this line:

       ```yaml
       password: password
       ```
       
       You can then use this password to log in to Code Server via your web browser.
 
