# How-to-access-NCBI-using-R
This tutorial will walk you through how to use R code to access NCBI and Pubmed
## Tutorial
### Step 1: Install and load the required packages
The following packages each have a specific function:
  ```rentrez``` is used to access the NCBI Entrez API for retrieving PubMed data.
  ```dplyr``` is used for data manipulation tasks.
  ```xml2``` is used to work with XML data.
```
install.packages("rentrez")
install.packages("dplyr")
install.packages("xml2")

library(rentrez)
library(dplyr)
library(xml2)
```

### Step 2: Retrive PubMed IDs (PMIDs)
To access PubMed articles, you typically start by retrieving a list of PubMed IDs (PMIDs) based on your search criteria. Here's an example of how to retrieve PMIDs using the ```entrez_search()``` function:
```
# Search for articles related to a specific query
search_query <- "your search query"

#Perform the search
search_results <- entrez_search(db = "pubmed", term = search_query, retmax = 5)

#Extract the IDs
pmids <- search_results$ids
print(pmids)
```
Replace "your search query" with your desired search terms, and retmax determines the maximum number of results to retrieve (in this case, 10). Print your results to view the PubMed IDs that were generated.
### Step 3: Fetch the MEDLINE records using IDs
Once you have the PMIDs, you can fetch the records using the entrez_fetch() function. Here's an example:
```
medline_records <- entrez_fetch(db, id = id_list, rettype = "medline", retmode = "xml")
```
### Step 4: Parse the XML data using 'xml2::read_xml'
```
parsed_records <- xml2::read_xml(medline_records)
```
### Step 5: Extract the titles, journals, and years from the parsed records
```
titles <- xml2::xml_text(xml2::xml_find_all(parsed_records, "//ArticleTitle"))
journals <- xml2::xml_text(xml2::xml_find_all(parsed_records, "//Journal"))
authors <- xml2::xml_text(xml2::xml_find_all(parsed_records, "//Author"))
years <- xml2::xml_text(xml2::xml_find_all(parsed_records, "//PubDate/Year"))
```
### Step 6: Extract the abstracts from the parsed records
```
abstracts <- xml2::xml_text(xml2::xml_find_all(parsed_records, "//AbstractText"))
```
### Step 7: Determine the Number of titles and subset the abstracts
This function will tell you how many titles were retrieved. This is good if you did not initally set the retmax() for a specific number.
```
num_titles <- length(titles)

# Subset the abstracts to match the number of titles
abstracts <- abstracts[1:num_titles]
```
### Step 8: Create the Article Info Data frame and display the results
```
article_info <- data.frame(
  Title = titles,
  Journal = journals,
  Year = years,
  Abstract = abstracts
)

print(article_info)
```

## More Data Exploration
### To count the occurrences of specific words
```
keyword1 <- "Your Keyword Here"
keyword2 <- "Your Keyword Here"

count_keyword1 <- sum(grepl(keyword1, abstracts, ignore.case = TRUE))
count_keyword2 <- sum(grepl(keyword2, abstracts, ignore.case = TRUE))

print(paste("Occurrences of", keyword1, "in abstracts:", count_keyword1))
print(paste("Occurrences of", keyword2, "in abstracts:", count_keyword2))
```
### Find the article with the longest title and print its details
The nchar function is used to calculate the number of characters in a string. By adding which.max, this will look for the longest titles only
```
longest_title_index <- which.max(nchar(titles))
longest_title <- titles[longest_title_index]
longest_title_authors <- authors[longest_title_index]
longest_title_year <- years[longest_title_index]
longest_title_abstract <- abstracts[longest_title_index]

print("Article with the longest title:")
print(paste("Title:", longest_title))
print(paste("Authors:", longest_title_authors))
print(paste("Year:", longest_title_year))
print(paste("Abstract:", longest_title_abstract))
```
### Determine the publication year with the highest number of articles and print the count:
```
year_counts <- table(years)
highest_year <- names(year_counts)[which.max(year_counts)]
highest_year_count <- max(year_counts)

print("Publication year with the highest number of articles:")
print(paste("Year:", highest_year))
print(paste("Count:", highest_year_count))
```
## Further Notes
XML data stands for eXtensible Markup Language much like HTML. It is designed to store and transport data.

In the code `xml2::read_xml`, `xml2::` is used to specify the `read_xml` function from the `xml2` package.

When multiple packages have functions with the same name, using `package_name::function_name` explicitly tells R which package's function to use. In this case, `read_xml` is a function from the `xml2` package, so we use `xml2::read_xml` to call that specific function.

By specifying `xml2::read_xml`, we ensure that we are using the `read_xml` function from the `xml2` package, even if there are other packages with functions named `read_xml`. This helps avoid any conflicts or ambiguity when multiple packages are loaded, as it explicitly tells R which package's function to use.

So, `xml2::read_xml` refers to the `read_xml` function from the `xml2` package, which is used to parse XML data in this code.
