/* login in to shake kemin

*/
desc tax_names;
+-------------+-------------+------+-----+---------+-------+
| Field       | Type        | Null | Key | Default | Extra |
+-------------+-------------+------+-----+---------+-------+
| id          | int(11)     |      | MUL | 0       |       |
| name        | text        | YES  | MUL | NULL    |       |
| unique_name | text        | YES  |     | NULL    |       |
| name_class  | varchar(20) | YES  |     | NULL    |       |
+-------------+-------------+------+-----+---------+-------+

select *
from tax_names
where name like 'Sporobolomyces roseus%';

+-------+-----------------------+-------------+-----------------+
| id    | name                  | unique_name | name_class      |
+-------+-----------------------+-------------+-----------------+
| 40563 | Sporobolomyces roseus |             | scientific name |
+-------+-----------------------+-------------+-----------------+


desc tax_nodes;
mysql:kemin:shake>desc tax_nodes;
+---------------------+-------------+------+-----+---------+-------+
| Field               | Type        | Null | Key | Default | Extra |
+---------------------+-------------+------+-----+---------+-------+
| id                  | int(11)     |      | PRI | 0       |       |
| parent              | int(11)     | YES  |     | NULL    |       |
| rank                | varchar(18) | YES  |     | NULL    |       |
| embl_code           | char(2)     | YES  |     | NULL    |       |
| division            | int(11)     | YES  |     | NULL    |       |
| inherited_div       | tinyint(1)  | YES  |     | NULL    |       |
| genetic_code        | int(11)     | YES  |     | NULL    |       |
| inherited_GC        | tinyint(1)  | YES  |     | NULL    |       |
| mtgenetic_code      | int(11)     | YES  |     | NULL    |       |
| inherited_MGC       | tinyint(1)  | YES  |     | NULL    |       |
| GenBank_hidden      | tinyint(1)  | YES  |     | NULL    |       |
| hidden_subtree_root | tinyint(1)  | YES  |     | NULL    |       |
| comments            | text        | YES  |     | NULL    |       |
+---------------------+-------------+------+-----+---------+-------+

select id, parent, rank, division, comments
from tax_nodes
where id=40563;
+-------+--------+---------+----------+----------+
| id    | parent | rank    | division | comments |
+-------+--------+---------+----------+----------+
| 40563 |   5429 | species |        4 |          |
+-------+--------+---------+----------+----------+


select *
from tax_names 
where id=5429;

select id, parent, rank, division, comments
from tax_nodes
where id=5429;

select id, parent, rank, division, comments

select id
from tax_nodes
where parent=5429;


