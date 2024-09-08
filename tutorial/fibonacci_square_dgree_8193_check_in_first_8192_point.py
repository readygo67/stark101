

from field import FieldElement
from fibonacci_square import *

'''
success generate trace
success build trace domain
success build polynomial
success build trace domain
success build lde domain!
success build lde evaluation!
success commit lde evaluation!
begin to build f_8193, cost too much time, be patient!
success build polynomial
f_8193.degree: 8192
check first 8192 points value are same!
Traceback (most recent call last):
  File "/Users/luokeep/Code/github.com/readygo67/stark101/tutorial/fibonacci_square_dgree_8193_check_in_first_8192_point.py", line 78, in <module>
    stark101_degree_8193()
  File "/Users/luokeep/Code/github.com/readygo67/stark101/tutorial/fibonacci_square_dgree_8193_check_in_first_8192_point.py", line 36, in stark101_degree_8193
    p0 = constraint_0(f_8193)
  File "/Users/luokeep/Code/github.com/readygo67/stark101/tutorial/fibonacci_square.py", line 89, in constraint_0
    p0 = numerator / 	denominator
  File "/Users/luokeep/Code/github.com/readygo67/stark101/tutorial/polynomial.py", line 209, in __truediv__
    assert mod == 0, 'Polynomials are not divisible.'
AssertionError: Polynomials are not divisible.
'''
def stark101_degree_8193_on_coset():
  trace_value = build_trace_eval()  # [y0,y1,y2,...,y1022]
  trace_domain = build_trace_domain()  # [g^0,g^1,g^2,...,g^1023]

  f = build_polynomial(trace_domain[:-1], trace_value)
  lde_domain = build_lde_domain(8192, True)
  lde_value = build_lde_value(f, lde_domain)
  lde_tree = commit(lde_value)
  assert lde_tree.root == '6c266a104eeaceae93c14ad799ce595ec8c2764359d7ad1b4b7c57a4da52be04'
  print('success commit lde evaluation!')

  lde_domain_8193 = lde_domain
  lde_domain_8193.append(lde_domain[-1]*FieldElement.generator()**2)
  lde_value_8193 = lde_value
  lde_value_8193.append(lde_value[-1] *2)
  assert len(lde_domain_8193) == 8193
  assert len(lde_value_8193) == 8193
  print('begin to build f_8193, cost too much time, be patient!')
  f_8193 = build_polynomial(lde_domain_8193, lde_value_8193)
  print('f_8193.degree:', f_8193.degree())

  for i in range (len(lde_domain_8193)-1):
    assert f.eval(lde_domain_8193[i]) == f_8193.eval(lde_domain_8193[i]) #check first 8192 points value
  print('check first 8192 points value are same!')


  channel = Channel()
  channel.send(lde_tree.root)

  p0 = constraint_0(f_8193)
  p1 = constraint_1(f_8193, trace_domain[-2])
  p2 = constraint_2(f_8193, trace_domain)

  alphas = channel.derive_alphas(3)
  cp = build_compostion_polynomial(p0, p1, p2, alphas[0], alphas[1], alphas[2])
  cp_eval = [cp(x) for x in lde_domain]

  cp_tree = commit(cp_eval)
  assert cp_tree.root == 'd7e5200e990727c6da6bf711aeb496244b8b48436bd6f29066e1ddb64e22605b'
  print('success commit cp evaluation!')

  channel.send(cp_tree.root)
  fri_polys, fri_domains, fri_layers, fri_merkles = fri_commit(cp, lde_domain, cp_eval,
                                                                    cp_tree, channel)
  print('success commit all proofs!')


  proof_0 = decommit_on_query(trace_domain, lde_domain, lde_value, lde_tree, fri_layers,
                              fri_merkles, 0, channel)
  valid = proof_0.verify(proof_0.final_values[0], lde_tree, fri_merkles)
  assert valid == True

  nb_pass = 0
  nb_fail = 0
  for i in range(8192-16):
   proof = decommit_on_query(trace_domain, lde_domain, lde_value, lde_tree, fri_layers,
                              fri_merkles, i, channel)
   valid = proof.verify(proof_0.final_values[0], lde_tree, fri_merkles)
   if valid == True:
    nb_pass += 1
   else:
    nb_fail += 1
    print(f'failed on index {i}')
  print(f'nb_pass: {nb_pass}, nb_fail: {nb_fail}')

'''
  使用8192插值太费时间，因此采用2048个点测试
'''
def stark101_test_2048_2048_on_coset():
  trace_value = build_trace_eval()  # [y0,y1,y2,...,y1022]
  trace_domain = build_trace_domain()  # [g^0,g^1,g^2,...,g^1023]

  f = build_polynomial(trace_domain[:-1], trace_value)
  lde_domain = build_lde_domain(2048, True)
  lde_value = build_lde_value(f, lde_domain)
  span = 2

  lde_tree = commit(lde_value)

  channel = Channel()
  channel.send(lde_tree.root)

  p0 = constraint_0(f)
  p1 = constraint_1(f, trace_domain[-2])
  p2 = constraint_2(f, trace_domain)

  alphas = channel.derive_alphas(3)
  cp = build_compostion_polynomial(p0, p1, p2, alphas[0], alphas[1], alphas[2])
  cp_eval = [cp(x) for x in lde_domain]

  cp_tree = commit(cp_eval)
  print(f'cp.degree: {cp.degree()}')

  channel.send(cp_tree.root)
  fri_polys, fri_domains, fri_layers, fri_merkles = fri_commit(cp, lde_domain, cp_eval, cp_tree,
                                                               channel)
  print('success commit all proofs!')

  proof_1 = decommit_on_query(trace_domain, lde_domain, lde_value, lde_tree, fri_layers,
                              fri_merkles, 1, channel, span)
  # proof_100 = decommit_on_query(trace_domain, lde_domain[:2048], lde_value, lde_tree, fri_layers,
  #                               fri_merkles, 100, channel, span)

  valid = proof_1.verify(proof_1.final_values[0], lde_tree, fri_merkles)
  assert valid == True
  # valid = proof_100.verify(proof_1.final_values[0], lde_tree, fri_merkles)
  # assert valid == True
  #
  # for i in range(10):
  #   idx = randint(0, 2048 - 2*span )
  #   proof = decommit_on_query(trace_domain, lde_domain[:2048], lde_value, lde_tree, fri_layers,
  #                             fri_merkles, idx, channel,span)
  #   valid = proof.verify(proof_1.final_values[0], lde_tree, fri_merkles)
  #   assert valid == True

  print(f'final values:', proof_1.final_values)
  print('success verify proofs')


'''
不在cost上做lde, 可能会造成verify 检查check cp0(x), f(x), f(g*x), f(g^2*x) 的关系时除0的错误
'''
def stark101_test_2048_2048_not_on_coset():
  trace_value = build_trace_eval()  # [y0,y1,y2,...,y1022]
  trace_domain = build_trace_domain()  # [g^0,g^1,g^2,...,g^1023]

  f = build_polynomial(trace_domain[:-1], trace_value)
  lde_domain = build_lde_domain(2048, False)
  lde_value = build_lde_value(f, lde_domain[:2048])
  span = 2
  for i in range(len(trace_value)):
    assert lde_domain[2*i] == trace_domain[i]
    assert lde_value[2*i] == trace_value[i]


  lde_tree = commit(lde_value)

  channel = Channel()
  channel.send(lde_tree.root)

  p0 = constraint_0(f)
  p1 = constraint_1(f, trace_domain[-2])
  p2 = constraint_2(f, trace_domain)

  alphas = channel.derive_alphas(3)
  cp = build_compostion_polynomial(p0, p1, p2, alphas[0], alphas[1], alphas[2])
  cp_eval = [cp(x) for x in lde_domain[:2048]]


  cp_tree = commit(cp_eval)

  channel.send(cp_tree.root)
  fri_polys, fri_domains, fri_layers, fri_merkles = fri_commit(cp, lde_domain[:2048], cp_eval, cp_tree,
                                                               channel)
  print('success commit all proofs!')

  proof_100 = decommit_on_query(trace_domain, lde_domain[:2048], lde_value, lde_tree, fri_layers,
                                fri_merkles, 100, channel, span)

  valid = proof_100.verify(proof_100.final_values[0], lde_tree, fri_merkles)
  assert valid == True

  for i in range(10):
    idx = randint(0, 2048 - 2*span )
    proof = decommit_on_query(trace_domain, lde_domain[:2048], lde_value, lde_tree, fri_layers,
                              fri_merkles, idx, channel,span)
    valid = proof.verify(proof_100.final_values[0], lde_tree, fri_merkles)
    assert valid == True

  print(f'final values:', proof_100.final_values)
  print('success verify proofs')


'''
扩展到4096,但是只用前2048个点，会造成trace_cp检查不成立，cp 各层之间的关系不成立。
'''
def stark101_test_2048_4096_on_coset():
  trace_value = build_trace_eval()  # [y0,y1,y2,...,y1022]
  trace_domain = build_trace_domain()  # [g^0,g^1,g^2,...,g^1023]

  f = build_polynomial(trace_domain[:-1], trace_value)
  lde_domain = build_lde_domain(4096, True)
  lde_value = build_lde_value(f, lde_domain[:2048])
  span = 2

  lde_tree = commit(lde_value)

  channel = Channel()
  channel.send(lde_tree.root)

  p0 = constraint_0(f)
  p1 = constraint_1(f, trace_domain[-2])
  p2 = constraint_2(f, trace_domain)

  alphas = channel.derive_alphas(3)
  cp = build_compostion_polynomial(p0, p1, p2, alphas[0], alphas[1], alphas[2])
  cp_eval = [cp(x) for x in lde_domain[:2048]]

  cp_tree = commit(cp_eval)

  channel.send(cp_tree.root)
  fri_polys, fri_domains, fri_layers, fri_merkles = fri_commit(cp, lde_domain[:2048], cp_eval, cp_tree,
                                                               channel)
  print('success commit all proofs!')

  proof_1 = decommit_on_query(trace_domain, lde_domain[:2048], lde_value, lde_tree, fri_layers,
                              fri_merkles, 1, channel, span)
  proof_100 = decommit_on_query(trace_domain, lde_domain[:2048], lde_value, lde_tree, fri_layers,
                                fri_merkles, 100, channel, span)

  valid = proof_1.verify(proof_1.final_values[0], lde_tree, fri_merkles)
  assert valid == True
  valid = proof_100.verify(proof_1.final_values[0], lde_tree, fri_merkles)
  assert valid == True

  for i in range(10):
    idx = randint(0, 2048 - 2*span )
    proof = decommit_on_query(trace_domain, lde_domain[:2048], lde_value, lde_tree, fri_layers,
                              fri_merkles, idx, channel,span)
    valid = proof.verify(proof_1.final_values[0], lde_tree, fri_merkles)
    assert valid == True

  print(f'final values:', proof_1.final_values)
  print('success verify proofs')



'''
  使用8192插值太费时间，因此采用2048个点测试
'''
def stark101_test_2048_2048_degree_on_coset():
  trace_value = build_trace_eval()  # [y0,y1,y2,...,y1022]
  trace_domain = build_trace_domain()  # [g^0,g^1,g^2,...,g^1023]

  f = build_polynomial(trace_domain[:-1], trace_value)
  lde_domain = build_lde_domain(2048, True)
  lde_value = build_lde_value(f, lde_domain)
  lde_tree = commit(lde_value)

  lde_domain_2050 = lde_domain
  lde_value_2050 = lde_value
  lde_domain_2050.append(FieldElement(1))   # fake f(g^0) = 1, to cheat verify
  lde_value_2050.append(FieldElement(1))
  lde_domain_2050.append(trace_domain[-2])   #fake =  f(g^1023) = 2338775057, to cheat verify
  lde_value_2050.append(trace_value[-1])
  lde_domain_2050.append(trace_domain[-3])   #fake =  f(g^1022) = xx, to cheat verify
  lde_value_2050.append(trace_value[-2])
  lde_domain_2050.append(trace_domain[-4])   #fake =  f(g^1021) = xx, to cheat verify
  lde_value_2050.append(trace_value[-3])
  lde_domain_2050.append(trace_domain[-2]*2)   #fake =  f(g^1023) = 2338775057, to cheat verify
  lde_value_2050.append(FieldElement(2338775057)*3)
  span = 2

  # assert len(lde_domain_2050) == 2050
  # assert len(lde_value_2050) == 2050
  print('begin to build f_2050, cost too much time, be patient!')
  f_2050 = build_polynomial(lde_domain_2050, lde_value_2050)
  print('f_2050.degree:', f_2050.degree())

  for i in range (len(lde_domain_2050)-2):
    assert f.eval(lde_domain_2050[i]) == f_2050.eval(lde_domain_2050[i]) #check first 8192 points value
  print('check first 2048 points value are same!')

  channel = Channel()
  channel.send(lde_tree.root)

  p0 = constraint_0(f_2050)
  p1 = constraint_1(f_2050, trace_domain[-2])
  p2 = constraint_2(f_2050, trace_domain)

  alphas = channel.derive_alphas(3)
  cp = build_compostion_polynomial(p0, p1, p2, alphas[0], alphas[1], alphas[2])
  cp_eval = [cp(x) for x in lde_domain]
  print('cp.degree:', cp.degree())

  cp_tree = commit(cp_eval)
  channel.send(cp_tree.root)
  fri_polys, fri_domains, fri_layers, fri_merkles = fri_commit(cp, lde_domain, cp_eval, cp_tree,
                                                               channel)
  print('success commit all proofs!')

  proof_1 = decommit_on_query(trace_domain, lde_domain, lde_value, lde_tree, fri_layers,
                              fri_merkles, 2, channel, span)
  proof_100 = decommit_on_query(trace_domain, lde_domain, lde_value, lde_tree, fri_layers,
                                fri_merkles, 100, channel, span)

  valid = proof_1.verify(proof_1.final_values[0], lde_tree, fri_merkles)
  assert valid == True
  valid = proof_100.verify(proof_1.final_values[0], lde_tree, fri_merkles)
  assert valid == True

  for i in range(10):
    idx = randint(0, 2048 - 2*span )
    proof = decommit_on_query(trace_domain, lde_domain, lde_value, lde_tree, fri_layers,
                              fri_merkles, idx, channel,span)
    valid = proof.verify(proof_1.final_values[0], lde_tree, fri_merkles)
    assert valid == True

  print(f'final values:', proof_1.final_values)
  print('success verify proofs')


'''
在lde_domain

'''
if __name__ == "__main__":
  # stark101_degree_8193_on_coset()
  #stark101_test_2048_2048_on_coset()
  stark101_test_2048_2048_degree_on_coset()
  #stark101_test_2048_2048_not_on_coset()

