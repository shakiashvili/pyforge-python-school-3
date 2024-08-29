import os
from dotenv import load_dotenv
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

# Load environment variables from .env file
load_dotenv()

# Read environment variables
DATABASE_USER = os.getenv("DB_USER", "user")
DATABASE_PASSWORD = os.getenv("DB_PASSWORD", "password")
DATABASE_NAME = os.getenv("DB_NAME", "db")
DATABASE_HOST = os.getenv("DB_HOST", "postgres")

# Construct the database URL
DATABASE_URL = os.getenv("DATABASE_URL", "postgresql://user:password@postgres:5432/db")

# Create the SQLAlchemy engine
engine = create_engine(DATABASE_URL)

# Create a configured "SessionLocal" class
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

# Create a base class for declarative models
Base = declarative_base()

# Dependency to get the database session


def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
