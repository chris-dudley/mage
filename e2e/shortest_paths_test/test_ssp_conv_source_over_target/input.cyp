MERGE (a:Node {id: 0}) MERGE (b:Node {id: 1}) MERGE (a)-[r:RELATION {id: 0}]->(b) ON CREATE SET r.weight = 2.0, r.capacity = 30.0;
MERGE (a:Node {id: 0}) MERGE (b:Node {id: 1}) MERGE (a)-[r:RELATION {id: 1}]->(b) ON CREATE SET r.weight = 1.0, r.capacity = 20.0;
MERGE (a:Node {id: 0}) MERGE (b:Node {id: 1}) MERGE (a)-[r:RELATION {id: 2}]->(b) ON CREATE SET r.weight = 3.0, r.capacity = 60.0;
MERGE (a:Node {id: 1}) MERGE (b:Node {id: 2}) MERGE (a)-[r:RELATION {id: 3}]->(b) ON CREATE SET r.weight = 2.0, r.capacity = 110.0;