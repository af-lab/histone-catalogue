# Installing Ubuntu on a USB stick
Instructions for installing Ubuntu with large persistent storage size on a bootable USB stick. It has been tested with Ubuntu 17.04 and 16.04.2 LTS and a 32 GB SanDisk USB 3.0 stick.

### mkusb
The installation uses mkusb which enables creation of a partition for persistent storage (casper-rw) to allow arbitrarily large space. Alternative tools such as unetbootin use a file for persistent storage and have a maximum of only 4 GB due to FAT limitations.

### Requirements
1. USB stick with 16-32 GB capacity to allow a persistant space of >8GB
1. Ubuntu iso file from https://www.ubuntu.com/download/desktop
1. Parent Ubuntu to run mkusb, which can be on another USB stick

## 1. Install mkusb on parent
Follow the instructions at https://help.ubuntu.com/community/mkusb

`$ sudo add-apt-repository ppa:mkusb/ppa`

`$ sudo apt update`

`$ sudo apt install mkusb mkusb-nox usb-pack-efi`

## 2. Create bootable Ubuntu on USB stick
1. Launch mkusb (eg click Dash and type `mkusb`)
1. Choose `i` = Install (make a boot device)
1. Choose `p` = Persistent live
1. Select the ubuntu iso file to use as source
1. Select the taget USB stick to create bootabel Ubuntu on
1. Check only `upefi` for persistent live settings
1. Pull slider for space for persistence choose value giving at least 10 GB (eg 50% of 32 GB USB stick). It is claimed that 100% can cause problems because there will be no space for swap files etc
1. Select `Go` and say OK to all the messages

## 3. Boot Ubuntu with USB stick
1. Adjust boot sequence on BIOS of PC to use USB as required
  * On Dell machines, hold F12 immediately on booting until "Preparing one-time boot menu" acknowledgement shows
1. Choose first UEFI USB device in list
1. Choose `Run Ubuntu (persistent live)` which is default
1. Booting can require a few minutes wait with a blank screen