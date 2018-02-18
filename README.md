# Tools



ADD DESCRIPTION






## NOTES for pushing and pulling linked packages

Prior to using package -- ensure that any changes on the linked side is pulled in. If changes were made to a linked package, please push the changes back up the subtree. 

Example: 
git subtree pull --prefix UMLS-Association https://github.com/henryst57/UMLS-Association.git master
git subtree push --prefix UMLS-Association https://github.com/henryst57/UMLS-Association.git master
