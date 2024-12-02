#!/usr/bin/env python3
""" simple rng implementation """

from Crypto.Cipher import AES
from Crypto.Hash import SHA256

class NIST_KAT_DRBG:
    """ AES-256 CTR to extract "fake" DRBG outputs that are compatible with
        the randombytes() call in the NIST KAT testing suite."""

    def __init__(self, seed):
        self.seed_length = 48
        assert len(seed) == self.seed_length
        self.key = b'\x00' * 32
        self.ctr = b'\x00' * 16
        update = self.get_bytes(self.seed_length)
        update = bytes(a^b for a,b in zip(update,seed))
        self.key = update[:32]
        self.ctr = update[32:]

    def __increment_ctr(self):
        x = int.from_bytes(self.ctr, 'big') + 1
        self.ctr = x.to_bytes(16, byteorder='big')

    def get_bytes(self, num_bytes):
        tmp = b''
        cipher = AES.new(self.key, AES.MODE_ECB)
        while len(tmp) < num_bytes:
            self.__increment_ctr()
            tmp  += cipher.encrypt(self.ctr)
        return tmp[:num_bytes]

    def random_bytes(self, num_bytes):
        output_bytes = self.get_bytes(num_bytes)
        update = self.get_bytes(48)
        self.key = update[:32]
        self.ctr = update[32:]
        return output_bytes
