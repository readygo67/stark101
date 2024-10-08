###############################################################################
# Copyright 2019 StarkWare Industries Ltd.                                    #
#                                                                             #
# Licensed under the Apache License, Version 2.0 (the "License").             #
# You may not use this file except in compliance with the License.            #
# You may obtain a copy of the License at                                     #
#                                                                             #
# https://www.starkware.co/open-source-license/                               #
#                                                                             #
# Unless required by applicable law or agreed to in writing,                  #
# software distributed under the License is distributed on an "AS IS" BASIS,  #
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.    #
# See the License for the specific language governing permissions             #
# and limitations under the License.                                          #
###############################################################################


from collections import Counter
from functools import reduce
from math import sqrt
from random import randint

from channel import Channel


def test_reproducability():
    c = Channel()
    c.send("Yes")
    r1 = c.receive_random_int(0, 2 ** 20)
    print(r1)
    d = Channel()
    d.send("Yes")
    r2 = d.receive_random_int(0, 2 ** 20)
    print(r2)
    assert r1==r2


def test_uniformity():
    range_size = 10
    c = Channel()
    c.send(str(randint(0, 2 ** 20)))
    num_tries = 2 ** 10
    dist = Counter(c.receive_random_int(0, range_size - 1) for i in range(num_tries))
    dist_rand = Counter(randint(0, range_size - 1) for i in range(num_tries))
    mean = num_tries / range_size
    normalized_stdv_channel = sqrt(sum((y-mean)**2 for y in dist.values())) / range_size
    normalized_stdv_random = sqrt(sum((y-mean)**2 for y in dist_rand.values())) / range_size
    assert abs(normalized_stdv_channel - normalized_stdv_random) < 4


if __name__ == "__main__":
    test_reproducability()
    test_uniformity()

