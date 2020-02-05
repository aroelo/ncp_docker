from project import db


class PavianInput(db.Model):
    """
CREATE TABLE pavian_data (
    file varchar(100),
    run varchar(100),
    sample varchar(100),
    nt varchar(10),
    date DATE,
    organism_taxid INT,
    organism_name varchar(100)
    support TINYINT(1)
    );
    """
    __tablename__ = "pavian_data"

    id = db.Column(db.Integer, primary_key=True)
    file = db.Column(db.String(100))
    run = db.Column(db.String(100))
    sample = db.Column(db.String(100))
    nt = db.Column(db.String(10))
    date = db.Column(db.Date)
    organism_taxid = db.Column(db.Integer)
    organism_name = db.Column(db.String(100))
    support = db.Column(db.Boolean)

    def __init__(self, file, run, sample, nt, date, organism_taxid, organism_name, support):
        self.file = file
        self.run = run
        self.sample = sample
        self.nt = nt
        self.date = date
        self.organism_taxid = organism_taxid
        self.organism_name = organism_name
        self.support = support