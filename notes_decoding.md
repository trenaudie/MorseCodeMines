# Notes and Learnings from coding the decoding script

### Opening a wav file 

The header takes 44 characters
file.read() takes as input a pointer to char, so we must use the address-of operator 
`&header[0]`. 
The raw file contents look like this 
```
RIFF�WAVEfmt D��Xdata�+Lv��"i�X�Y��L��Bj����[(��5�I�%�L�}�������%�B�b�������#�l���*��+���~�H�*�%�:�i����.�������������h���<���A�e��(  �
=
$�����W�.bnS�
RIFF�WAVEfmt D��Xdata�+Lv��"i�X�Y��L��Bj����[(��5�I�%�L�}�������%�B�b�������#�l���*��+���~�H�*�%�:�i����.�������������h���<���A�e��( �
=
$�����W�.bnS�
```
This is obviously not encoded as a string, but as binary data (01101000 00011011 etc ), where each pair of 2 bytes corresponds to an audio sample.
We must convert these 2 byte pairs into their float values, hence the code 
```cpp
short value = *(short*)buffer;
float normalized = value / 32768.0f;
```


