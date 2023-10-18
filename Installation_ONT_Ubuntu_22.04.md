#This is a file created out of frustration, while trying to install ONT in ubuntu 22.04

We are following the instructions placed in the [ONT_how_to_repository](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/experiment-companion-minknow/v/mke_1013_v1_revcy_11apr2016/installing-minknow-on-linu)

##For Ubuntu 20:

```bash
sudo apt update
sudo apt install wget
wget -O- https://cdn.oxfordnanoportal.com/apt/ont-repo.pub | sudo apt-key add -
echo "deb http://cdn.oxfordnanoportal.com/apt focal-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
```

This code breaks apart because apt-key has been deprecated in ubuntu 22. Instead we need to use gpg and follow a different protocol.

For the modified protocol update the system as before (it is just always a good idea to do so), and install wget if you don't have done so (I doubt it, but check it exists first)

```bash
sudo apt update
sudo apt install wget
```

Now the fun begins. Download the public key from ONT and save it as a file (omit the -O- flag from wget)

```bash
wget https://cdn.oxfordnanoportal.com/apt/ont-repo.pub
```

Following that we need to get the key in a readable format for gpg. Verify that the filetype is "PGP public key block Public-Key (old)". You are going to type:

```bash
file ont-repo.pub
```

and you should get the "ont-repo.pub: PGP public key block Public-Key (old)" output. If that is the case proceed to the following step. 
Now, if you start lookinga the new way keys will be managed in ubuntu, you will realize that gpg supports a number of key formats, so if your key is in a different format, convert it by importing it into a temp keyring (prior to continuing to the next step), then exporting it again:

```bash
gpg --no-default-keyring --keyring ./temp-keyring.gpg --import ont-repo.pub
gpg --no-default-keyring --keyring ./temp-keyring.gpg --export --output ont-repo.gpg
rm temp-keyring.gpg
```

if you are unceratin if your key is in the right format, just run the previous command and add the appropriate extension to it. Also, if a ```rm temp-keyring.gpg~``` file is created, you can delete it as well.

After you convert the key, **do not** add it to apt's trusted keystore (/etc/apt/trusted.gpg.d/). what you should do now is to put it in /etc/apt/keyrings/. Create that keyrings directory if it doesn't exist, first. That location is a convention recommended by ```man 5 sources.list``` in Ubuntu 22.04.

```bash
sudo mv ont-repo.gpg /etc/apt/keyrings/
```
Now, I borrow from the community knowledge (see more [here](https://askubuntu.com/questions/1286545/what-commands-exactly-should-replace-the-deprecated-apt-key)). According to the previous source: "At this point, nothing has changed and ```apt``` doesn't know the key exists. The last step is to modify the specific ```.list``` file for the repository to tell apt where to find the key for that specific repo."

Edit the file ```/etc/apt/sources.list.d/nanoporetech.sources.list```, and in between ```deb``` and the url, add ```[signed-by=/etc/apt/keyrings/ont-repo.gpg]```. Use your favorite editor and remember that this file can only be edited with sudo permissions (very likely).

My favorite way to do this will be:

```bash
echo "deb [signed-by=/etc/apt/keyrings/ont-repo.gpg] http://cdn.oxfordnanoportal.com/apt focal-stable non-free" | sudo tee /etc/apt/sources.list.d/nanoporetech.sources.list
```

This will ensure that apt recognizes the ont repo for installation...  now the test of fire! 

Now we can go back to our ONT friendly support group and continue with the installation

##Install MinKNOW using the command:
###CPU version:

```bash
sudo apt update
sudo apt install ont-standalone-minknow-release
```

###GPU version:

```bash
sudo apt update
sudo apt install ont-standalone-minknow-gpu-release
```

For GPU specs and a list of example GPU models, refer to the [MinION Mk1B IT requirements document](https://community.nanoporetech.com/requirements_documents/minion-it-reqs.pdf).

Now you can go back to the [ONT support page](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/experiment-companion-minknow/v/mke_1013_v1_revcy_11apr2016/installing-minknow-on-linu) and continue checking that the installation is correct and everything is in order

