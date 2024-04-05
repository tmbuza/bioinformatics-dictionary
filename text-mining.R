# Load the tm package
library(tm)

# Sample text for analysis
text <- c(
"Programming experience with two or more programming languages including: Python, R for bioinformatics data analysis.
Proficiency in working with bulk and single cells NGS data, bioinformatics databases including PPI, CRISPR, Proteomics, or Metabolome.
Experience in biological pathways analysis, systems network biology a plus
Experience in machine learning, especially deep learning a plus
Two years hand on experience of bioinformatics data analysis.
Good understand of Cancer Biology, Metabolisms, and Immunology.
Ability to develop and benchmark machine learning Client algorithms.
Experience in software and/or pipeline development a plus."
)

# Create a Corpus object from the text
corpus <- Corpus(VectorSource(text))

# Preprocessing: Convert to lowercase, remove punctuation, and remove stopwords
corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeWords, stopwords("en"))

# Stemming: Reduce words to their root form
corpus <- tm_map(corpus, stemDocument)

# Create a TermDocumentMatrix to represent the document-term matrix
dtm <- DocumentTermMatrix(corpus)

# Convert the document-term matrix to a matrix
matrix <- as.matrix(dtm)

# Calculate the frequency of terms
term_freq <- colSums(matrix)

# Display the top 5 most frequent terms
top_terms <- head(sort(term_freq, decreasing = TRUE), 5)
print("Top 5 most frequent terms:")
print(top_terms)
