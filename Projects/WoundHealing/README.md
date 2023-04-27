# Wound Healing

April 10, 2023

---

- GL: Great Lakes
- LS: Lab Server
- Lab: Wet Lab

## DATA PIPELINE
1. (Lab) Jillian and Walter set experiments
2. (Lab) Experiment documentation written
3. (Lab) Zeiss generates images over time
4. (LS) Images are moved to lab server
5. (Lab) Documentation of experiment completion
6. (GL) Images are moved to Great Lakes

## PROCESSING PIPELINE
1. (GL) Image Preprocessing
2. (GL) Cell Tracking
3. (GL) Hypergraph Analysis

## Meeting Notes

### April 20, 2023
- **Indika:** He thinks we don't need Hooks
- **Cooper:** Stardist does a great job at cell segmentation (Cell segmentation), he thinks we need Hooks to identify the nuclei, count the effectivness of the PIP-FUCCI assay, etc.
- **Jillian:** (data generation) New set of experiments, imaging and stitching on Zeiss, put on turbo
- **Joshua:** Cell tracking (from Cooper on turbo: tensors that are time x image height x image width with *i,j,k* coordinate indicating if a cell is present)

### April 17, 2023
Progress:
1. Amit has been sent only 3 pictures

TODO:
1. Set up pipeline for small amounts of data to be shared
2. Set up pipeline to get Amit's code and run it on GreatLakes
3. 

### April 10, 2023

Problems:
1. Monitor wound boundary over time
2. Monitor individual cell migration

Pipeline:
1. (lab) Jillian dos wicked cool science in the lab
2. (lab) The microscope runs for 3-4 days
3. (zeiss) Images are stiched together on Zeiss (Output: scenes x channels x time)
4. () Preprocess images


Timeline:
- Today: Oblique, PIP-FUCCI
- 
