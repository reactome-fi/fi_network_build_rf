This project is the updated version of the FI network construction project. The original project, which is hosted at [Reactome's FI Network Project](https://github.com/reactome/fi_network_build), is based on training a Naive Bayes Classifier. This updated version is based on training random forest. There are two functions of this project:

* Generate FI files from Reactome and other manually curated pathway databases so that they can be used to generated training and test matrix files to train random forest. The project uses these files is [fi-network-ml at reactome-idg](https://github.com/reactome-idg/fi-network-ml). 


* Integrate the predicted results from the fi-network-ml project together with extracted pathway FIs to develop a FI network database and Cytoscape files so that they can be used to update ReactomeFIViz.

Note: Some of jar files need to be installed locally. Use ant scripts in the ant folder and jar files in the install_jar to do that. Please note the license statement for jacksum.jar.