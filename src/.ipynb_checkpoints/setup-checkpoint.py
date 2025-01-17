from setuptools import setup, find_packages

setup(
    name="verkkoFillet",                 # Your package name
    version="0.1.0",                      # Version
    author="Juhyun Kim",                   # Author name
    author_email="kimj75@nih.gov", # Author email
    description="A toolkit for cleaning Verkko assemblies.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jjuhyunkim/verkko-fillet", # Repository URL
    packages=find_packages(),             # Automatically find packages in your directory
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[                    # Dependencies
        "numpy", 
        "pandas",
        
        
    ],
)
