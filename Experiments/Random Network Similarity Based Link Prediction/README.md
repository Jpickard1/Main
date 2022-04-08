# Random Network Similarity Based Link Prediction

Auth: Joshua Pickard

Date: April 8, 2022

---

This experiment tests the effectiveness of various similarity indices used for link prediction on random networks. Insipred by the papers "Finding missing edges and communities in incomplete networks" by Yang and Gregory, which poses 4 reasons edges may be missing from a network, and the paper "Uncovering missing links with cold ends" by Zhu, which suggests LHN is an important similarity index for link prediction when a network has certain characteristics, this experiment tests 13 siumilarity indices (10 local and 3 global) for link prediction on different types of networks.

Four random types of networks are generated

Four different schemes of edge removal are applied

I want to look at how the eigenspectrum changes for each type of graph under each type of edge removal.
