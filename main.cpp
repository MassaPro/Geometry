#include <iostream>

using namespace std;

struct Node {
  Node *next = nullptr;
  int data = 0;
};

struct Queue {
  Node *beg = new Node();
  Node *end = new Node();
  size_t size = 0;
};

void push(Queue &q, int data) {
  if (!q.size) {
    Node *new_node = new Node;
    new_node->data = data;

    q.end = new_node;
    q.beg = new_node;
  } else {
    Node *new_end = new Node;
    new_end->data = data;
    q.end->next = new_end;
    q.end = new_end;
  }

  q.size++;
}

int pop(Queue &q) {
  int beg_old = 0;
  if (q.size > 0) {
    beg_old = q.beg->data;

    Node *old_beg = q.beg;
    q.beg = q.beg->next;

    delete old_beg;
    q.size--;
  }

  return beg_old;
}

void priv_push(Queue &q, int data) {
  if (q.size > 1) {
    Node *new_node = new Node;
    new_node->data = data;

    Node *where_node = q.beg;

    size_t site = (q.size + 1) / 2;

    for (size_t i = 0; i < site; i++) {
      where_node = where_node->next;
    }
    new_node->next = where_node->next;

    where_node->next = new_node;

    q.size++;
  } else {
    push(q, data);
  }
}

void clear(Queue &q) {
  while (q.size) {
    pop(q);
  }
}

int main() {
  Queue queue_goblins;
  int n;
  cin >> n;

  char elem;
  int number;

  for (int i = 0; i < n; i++) {
    cin >> elem;

    if (elem == '+') {
      cin >> number;
      cout << i << " " << number;
      push(queue_goblins, number);
    }
    if (elem == '*') {
      cin >> number;
      cout << i << " " << number;
      priv_push(queue_goblins, number);
    }
    if (elem == '-') {
      cout << i << " ";
      cout << pop(queue_goblins) << '\n';
    }
  }

  clear(queue_goblins);

  return 0;
}