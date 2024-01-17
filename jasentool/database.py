"""Module for handling mongodb requests"""
import pymongo

class Database:
    """Class that assists in handling mongodb request"""
    uri = "mongodb://localhost:27017/"
    db = None

    @staticmethod
    def initialize(db_name):
        """Initialize mongodb client"""
        client = pymongo.MongoClient(Database.uri)
        Database.db = client[db_name] # Database Name
        Database.db_name = db_name # Database Name
        #Database.collection = Database.db["sample"] # Collection Name

    @staticmethod
    def insert(collection, data):
        """Insert data into mongodb"""
        Database.db[collection].insert(data)

    @staticmethod
    def find(collection, query, fields):
        """Find data in mongodb"""
        return Database.db[collection].find(query, fields)

    @staticmethod
    def find_one(collection, query):
        """Find one entry in mongodb"""
        return Database.db[collection].find_one(query)

    @staticmethod
    def get_pvl(collection, query):
        """Get pvl result data from mongodb"""
        return Database.db[collection].find(query, {"_id": 0, "aribavir.lukS_PV.present": 1})

    @staticmethod
    def get_mlst(collection, query):
        """Get mlst result data from mongodb"""
        return Database.db[collection].find(query, {"_id": 0, "mlst": 1})

    @staticmethod
    def get_cgmlst(collection, query):
        """Get cgmlst result data from mongodb"""
        return Database.db[collection].find(query, {"_id": 0, "alleles": 1})

    @staticmethod
    def get_meta_fields():
        """Get respective metadata from mongodb"""
        fields = {
            "id": 1,
            "mlst.sequence_type": 1,
            "aribavir.lukF_PV.present": 1,
            "aribavir.lukS_PV.present": 1,
            "missing": 1,
            "metadata.QC": 1,
            "metadata.Comment": 1,
            "run": 1
        }
        return fields
