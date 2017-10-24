SELECT  `document`.`id`,
        `document`.`title`,
        COUNT(`comment`.`id`) as `comments`

FROM `document_document` as `document`

LEFT JOIN `django_comments` as `comment`
  ON `comment`.`object_pk` = `document`.`id`
    AND `comment`.`content_type_id` = (SELECT `id` FROM `django_content_type` WHERE `django_content_type`.`model` = 'document')

GROUP BY `document`.`id`

ORDER BY `comments` DESC

LIMIT 100
