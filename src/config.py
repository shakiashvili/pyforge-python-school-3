# import os
# from pydantic_settings import BaseSettings


# class Settings(BaseSettings):
#     DATABASE_URL: str = os.getenv
#     ("DATABASE_URL",
#      "postgresql+psycopg2://postgres:password@localhost/mydatabase")
#     ENVIRONMENT: str = os.getenv("ENVIRONMENT", "development")


# settings = Settings()
# print(settings.DATABASE_URL)

from pydantic_settings import BaseSettings


class Settings(BaseSettings):
    database_url: str  # Make sure this matches your environment variable name

    class Config:
        env_file = ".env"  # Path to your .env file if you have one
        env_file_encoding = 'utf-8'  # Ensure proper encoding if necessary

# Create a settings instance


settings = Settings()
