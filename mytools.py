

#partition a list into sublists of len n
def partn(l,n):
  nl = len(l)
  if nl%n == 0:
    nparts =  nl/n
  else:
    nparts = nl/n + 1
  return [l[i*n:(i+1)*n] for i in range(nparts)]


#round-up division of integers
def mydiv(n,p):
  if n%p == 0:
    return n/p
  else:
    return n/p+1

#linked list
class linked_node:
  def __init__(self,data):
    self.prev = None
    self.next = None
    self.data = data

class linked_list:
  def __init__(self):
    self.curr = None
    self.head = None
    self.tail = None
  
  def push_back(self,data):
    newtail = linked_node(data)
    if self.tail is not None:
      newtail.prev = self.tail
      self.tail.next = newtail
      self.tail = newtail
    else:
      self.head = newtail
      self.tail = newtail
      self.curr = newtail

  def insert(self,data):
    #inserts a node after the current node
    #returns the inserted node
    newnode = linked_node(data)
    if self.curr is not None:
      if self.curr.next is not None:
        self.curr.next.prev = newnode
        newnode.next = self.curr.next
        newnode.prev = self.curr
        self.curr.next = newnode
        self.curr = newnode
      else: #adding to tail of list
        self.curr.next = newnode
        newnode.prev = self.curr
        self.curr = newnode
        self.tail = newnode
    else:
      self.head = newnode
      self.tail = newnode
      self.curr = newnode
    return self.curr

  def delete(self):
    #deletes the current node
    #returns the next node
    if self.curr is not None:
      if self.curr.prev is not None:
        self.curr.prev.next = self.curr.next
      else: #deleting the head
        self.head = self.curr.next
      if self.curr.next is not None:
        self.curr.next.prev = self.curr.prev
      else: #deleting the tail
        self.tail = self.curr.prev
      self.curr = self.curr.next
    return self.curr

  def step(self):
    self.curr = self.curr.next

  def tohead(self):
    self.curr = self.curr.head

  def printll(self):
    self.curr = self.head
    while self.curr is not None:
      print self.curr.data
      self.step()
    self.curr = self.head
