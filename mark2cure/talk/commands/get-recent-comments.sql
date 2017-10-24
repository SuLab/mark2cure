SELECT  `comment`.`user_id`,
        `comment`.`user_name`,
        `comment`.`comment`,
        `comment`.`submit_date`,
        `document`.`document_id` as `pmid`

FROM `django_comments` as `comment`

INNER JOIN `document_document` as `document`
  ON `document`.`id` = `comment`.`object_pk`

WHERE `comment`.`content_type_id` = (SELECT `id` FROM `django_content_type` WHERE `django_content_type`.`model` = 'document')

ORDER BY `comment`.`submit_date` DESC

LIMIT 100
